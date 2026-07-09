// Copyright 2021-2026 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::process;
use std::time::Instant;

use anyhow::{anyhow, Context};
use clap::Parser;
use cli::Cli;
use log::{error, info, warn};

mod cli;
mod demux;
mod utils;

fn main() -> anyhow::Result<()> {
    let start_time = Instant::now();
    let cli = Cli::parse();

    utils::setup_logging(cli.quiet)?; // Settting up logging
    let forward_format = utils::which_format(&cli.forward)?;
    let mut output_format = forward_format;
    let mismatch = cli.mismatch;
    let level = utils::to_niffler_level(cli.level);

    // Handle compression
    if let Some(fmt_str) = cli.format {
        output_format = utils::to_niffler_format(fmt_str)?;
        info!(
            "Output files will {} compressed",
            utils::to_compression_ext(output_format)
        );
    }

    // Determine mode
    let is_pe = cli.reverse.is_some();
    info!(
        "sabreur v0.7 starting up in {} mode",
        if is_pe { "paired-end" } else { "single-end" }
    );

    // Read and validate the barcode file up front, before touching the
    // output directory, so a malformed barcode file fails fast without
    // side effects on disk.
    let barcode_content = fs::read_to_string(&cli.barcode)?;
    let barcode_fields = utils::split_by_tab(&barcode_content)?;
    utils::validate_barcode_fields(&barcode_fields, is_pe)
        .with_context(|| format!("Invalid barcode file '{}'", cli.barcode))?;
    let mut barcode_info: demux::Barcode = HashMap::new();
    let mut record_stats: HashMap<&[u8], u32> = HashMap::new();

    if mismatch != 0 {
        warn!("Allowing up to {} mismatches", mismatch);
    }

    // Output directory handling
    if cli.output.exists() {
        if !cli.force {
            error!(
                "Output folder '{}' already exists. Use --force to overwrite.",
                cli.output.display()
            );
            process::exit(exitcode::CANTCREAT);
        }

        info!("Reusing directory {}", cli.output.display());
        fs::remove_dir_all(&cli.output).with_context(|| {
            anyhow!(
                "Failed to remove '{}'. Check your permissions.",
                cli.output.display()
            )
        })?;
    }

    fs::create_dir_all(&cli.output).with_context(|| {
        anyhow!(
            "Failed to create '{}'. Check your permissions.",
            cli.output.display()
        )
    })?;

    // Helper to create a persistent, per-output-file writer. Each writer is
    // opened and wrapped with its compressor exactly once here, then reused
    // for every record written to that file, this is what avoids
    // re-creating a compressor per record. Return type must match
    // utils::create_output_writer exactly: Box<dyn Write + Send>, wrapped
    // in anyhow::Result (not niffler::Error -- anyhow converts it via `?`
    // inside create_output_writer already).
    let create_writer =
        |name: &str, format| -> anyhow::Result<Box<dyn std::io::Write + Send + 'static>> {
            let path = utils::create_relpath_from(&cli.output, name, format);
            utils::create_output_writer(&path, format, level)
        };

    // Main processing
    if let Some(reverse_path) = &cli.reverse {
        let mut reverse_format = utils::which_format(reverse_path)?;
        if output_format != niffler::send::compression::Format::No {
            reverse_format = output_format;
        }

        for fields in &barcode_fields {
            let forward_writer = create_writer(fields[1], output_format)?;
            let reverse_writer = create_writer(fields[2], reverse_format)?;
            barcode_info.insert(fields[0].as_bytes(), vec![forward_writer, reverse_writer]);
        }

        let unknown_fwd_path =
            utils::create_relpath_from(&cli.output, "unknown_R1.fa", output_format);
        let unknown_rev_path =
            utils::create_relpath_from(&cli.output, "unknown_R2.fa", reverse_format);
        let unknown_fwd = utils::create_output_writer(&unknown_fwd_path, output_format, level)?;
        let unknown_rev = utils::create_output_writer(&unknown_rev_path, reverse_format, level)?;
        barcode_info.insert(b"XXX", vec![unknown_fwd, unknown_rev]);

        let (stats, unk_status) = demux::pe_demux(
            &cli.forward,
            reverse_path,
            &mut barcode_info,
            mismatch,
            &mut record_stats,
        )?;

        if !cli.quiet {
            for (barcode, count) in stats {
                info!(
                    "{} records found for barcode {}",
                    count,
                    String::from_utf8_lossy(barcode)
                );
            }
        }

        // Flush every writer now, before deciding which (if any) unknown
        // files to delete. Writers are held open for the whole run for
        // performance, so nothing is guaranteed to be fully on disk until
        // this point.
        for files in barcode_info.values_mut() {
            for writer in files.iter_mut() {
                writer.flush().context("Failed to flush an output file")?;
            }
        }

        match unk_status.as_str() {
            "truetrue" => {
                fs::remove_file(&unknown_fwd_path)?;
                fs::remove_file(&unknown_rev_path)?;
            }
            "truefalse" => fs::remove_file(&unknown_fwd_path)?,
            "falsetrue" => fs::remove_file(&unknown_rev_path)?,
            _ => {}
        }
    } else {
        for fields in &barcode_fields {
            let writer = create_writer(fields[1], output_format)?;
            barcode_info.insert(fields[0].as_bytes(), vec![writer]);
        }

        let unknown_path = utils::create_relpath_from(&cli.output, "unknown.fa", output_format);
        let future_unk = unknown_path.clone();
        let unknown_writer = utils::create_output_writer(&unknown_path, output_format, level)?;
        barcode_info.insert(b"XXX", vec![unknown_writer]);

        let (stats, unk_empty) =
            demux::se_demux(&cli.forward, &mut barcode_info, mismatch, &mut record_stats)?;

        if !cli.quiet {
            for (barcode, count) in stats {
                info!(
                    "{} records found for barcode {}",
                    count,
                    String::from_utf8_lossy(barcode)
                );
            }
        }

        // Flush every writer now that all records have been written. See
        // the equivalent comment in the paired-end branch above.
        for files in barcode_info.values_mut() {
            for writer in files.iter_mut() {
                writer.flush().context("Failed to flush an output file")?;
            }
        }

        if unk_empty {
            fs::remove_file(future_unk)?;
        }
    }

    if !cli.quiet {
        let duration = start_time.elapsed();
        info!("Results saved in {}", cli.output.display());
        info!(
            "Walltime: {}h:{}m:{}s {}ms",
            duration.as_secs() / 3600,
            (duration.as_secs() / 60) % 60,
            duration.as_secs() % 60,
            duration.as_millis() % 1000
        );
        info!("Thanks. Share. Come again!");
    }

    Ok(())
}
