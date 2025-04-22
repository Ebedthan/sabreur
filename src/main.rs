// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::{self, OpenOptions};
use std::process;
use std::time::Instant;

use anyhow::{anyhow, Context};
use clap::Parser;
use cli::Cli;
use log::{error, info, warn};

mod cli;
mod demux;
mod utils;

// TODO: Check if supplied barcode file for se or pe is properly
// formated before giving it to the demultiplexing function
fn main() -> anyhow::Result<()> {
    let start_time = Instant::now();
    let cli = Cli::parse();

    utils::setup_logging(cli.quiet)?; // Settting up logging
    let forward_format = utils::which_format(&cli.forward);
    let mut output_format = forward_format;
    let mismatch = cli.mismatch;
    let raw_level = cli.level;

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

    // Read barcode file
    let barcode_content = fs::read_to_string(&cli.barcode)?;
    let barcode_fields = utils::split_by_tab(&barcode_content)?;
    let mut barcode_info: demux::Barcode = HashMap::new();
    let mut record_stats: HashMap<&[u8], u32> = HashMap::new();

    if mismatch != 0 {
        warn!("Allowing up to {} mismatches", mismatch);
    }

    // Helper to create writer
    let create_writer = |name: &str, format| -> anyhow::Result<_> {
        let path = utils::create_relpath_from(&cli.output, name, format);
        Ok(OpenOptions::new().create(true).append(true).open(path)?)
    };

    // Main processing
    if let Some(reverse_path) = &cli.reverse {
        let mut reverse_format = utils::which_format(reverse_path);
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
        let unknown_fwd = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&unknown_fwd_path)?;
        let unknown_rev = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&unknown_rev_path)?;
        barcode_info.insert(b"XXX", vec![unknown_fwd, unknown_rev]);

        let (stats, unk_status) = demux::pe_demux(
            &cli.forward,
            reverse_path,
            output_format,
            utils::to_niffler_level(raw_level),
            &barcode_info,
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
        let unknown_writer = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&unknown_path)?;
        barcode_info.insert(b"XXX", vec![unknown_writer]);

        let (stats, unk_empty) = demux::se_demux(
            &cli.forward,
            output_format,
            utils::to_niffler_level(raw_level),
            &barcode_info,
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
