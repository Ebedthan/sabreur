// Copyright 2021-2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate anyhow;
extern crate clap;
extern crate exitcode;
extern crate sysinfo;

use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::process;
use std::time::Instant;

use anyhow::{anyhow, Context};
use clap::crate_version;
use log::{error, info, warn};
use sysinfo::System;

mod app;
mod demux;
mod utils;

// TODO: Check if supplied barcode file for se or pe is properly
// formated before giving it to the demultiplexing function
fn main() -> anyhow::Result<()> {
    let mut sys = System::new_all();
    sys.refresh_all();
    let startime = Instant::now();

    // Define command-line arguments ----------------------------------------
    let matches = app::build_app().get_matches_from(env::args_os());

    // is --quiet option specified by the user?
    let quiet = matches.get_flag("quiet");
    utils::setup_logging(quiet)?; // Settting up logging

    // Read command-line arguments
    let forward = matches
        .get_one::<String>("FORWARD")
        .expect("input file is required");

    let mut forward_format = utils::which_format(forward);

    let barcode = matches
        .get_one::<String>("BARCODE")
        .expect("input barcode is required");

    let output: &PathBuf = matches.get_one("output").unwrap();
    let mismatch: u8 = *matches.get_one("mismatch").unwrap();

    // If user force output to be compressed even if input is not
    // add option to change compression of output
    let mut format = niffler::send::compression::Format::No;
    if matches.contains_id("format") {
        format = utils::to_niffler_format(
            matches.get_one::<String>("format").unwrap(),
        )
        .with_context(|| {
            anyhow!("Could not convert compression format to niffler format")
        })?;
    }

    let raw_level: u8 = *matches.get_one("level").unwrap();
    let force = matches.get_flag("force");

    info!("sabreur v{} starting up!", crate_version!());
    if matches.contains_id("REVERSE") {
        info!("You are in single-end mode");
    } else {
        info!("You are in paired-end mode");
    }

    // Change file compression format here for files extension
    if format != niffler::send::compression::Format::No {
        forward_format = format;
        info!(
            "Output files will be {} compressed",
            utils::to_compression_ext(forward_format)
        );
    }

    // Handle output dir
    let outdir_exists = output.exists();
    if outdir_exists && !force {
        error!("output folder '{}', already exists! change it using --out or use --force", output.display());
        process::exit(exitcode::CANTCREAT);
    } else if outdir_exists && force {
        info!("Reusing directory {}", output.display());
        fs::remove_dir_all(output).with_context(|| anyhow!("Could not remove folder '{}'. Do you have permission to remove this folder?", output.display()))?;
        fs::create_dir(output).with_context(|| anyhow!("Could not create folder '{}'. Do you have permission to create this folder?", output.display()))?;
    } else if !outdir_exists {
        fs::create_dir(output)?;
    }

    // Read data from barcode file
    let mut barcode_info: demux::Barcode = HashMap::new();
    let barcode_data = fs::read_to_string(barcode)?;
    let barcode_fields = utils::split_by_tab(&barcode_data).unwrap();

    if mismatch != 0 {
        warn!("Barcode mismatch allowed: {}", mismatch);
    }

    let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

    // Main processing of reads
    match !matches.contains_id("REVERSE") {
        // single-end fasta mode
        true => {
            let ext = utils::to_compression_ext(forward_format);
            // Read barcode data
            for b_vec in barcode_fields.iter() {
                let file_path = utils::create_relpath_from(
                    &mut output.clone(),
                    [format!("{}{}", b_vec[1], ext).as_str()].to_vec(),
                );
                let file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(file_path)?;
                barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
            }
            // Create unknown file
            let unk_path = utils::create_relpath_from(
                &mut output.clone(),
                [format!("{}{}", "unknown.fa", ext).as_str()].to_vec(),
            );
            let future_unk_path = unk_path.clone();
            let unknown_file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(unk_path)?;
            barcode_info.insert(b"XXX", vec![unknown_file]);
            // Demultiplexing
            let (stats, is_unk_empty) = demux::se_demux(
                forward,
                format,
                utils::to_niffler_level(raw_level),
                &barcode_info,
                mismatch,
                &mut nb_records,
            )?;
            if !quiet {
                for (key, value) in stats.iter() {
                    info!(
                        "{} records found for {} barcode",
                        value,
                        String::from_utf8_lossy(key)
                    );
                }
            }
            if is_unk_empty {
                fs::remove_file(future_unk_path)?;
            }
        }
        // paired-end fasta mode
        false => {
            let reverse = matches.get_one::<String>("REVERSE").unwrap();
            let mut reverse_format = utils::which_format(reverse);
            if format != niffler::send::compression::Format::No {
                reverse_format = format;
            }
            let f_ext = utils::to_compression_ext(forward_format);
            let r_ext = utils::to_compression_ext(reverse_format);
            // Read barcode data
            for b_vec in barcode_fields.iter() {
                let file_path1 = utils::create_relpath_from(
                    &mut output.clone(),
                    [format!("{}{}", b_vec[1], f_ext).as_str()].to_vec(),
                );
                let file_path2 = utils::create_relpath_from(
                    &mut output.clone(),
                    [format!("{}{}", b_vec[2], r_ext).as_str()].to_vec(),
                );
                let file1 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(file_path1)?;
                let file2 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(file_path2)?;
                barcode_info.insert(b_vec[0].as_bytes(), vec![file1, file2]);
            }
            // Create unknown files
            let unk_path1 = utils::create_relpath_from(
                &mut output.clone(),
                [format!("{}{}", "unknown_R1.fa", f_ext).as_str()].to_vec(),
            );
            let unk_path2 = utils::create_relpath_from(
                &mut output.clone(),
                [format!("{}{}", "unknown_R2.fa", r_ext).as_str()].to_vec(),
            );
            let future_unk_path1 = unk_path1.clone();
            let future_unk_path2 = unk_path2.clone();
            let unknown_file1 = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(unk_path1)?;
            let unknown_file2 = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(unk_path2)?;
            barcode_info.insert(b"XXX", vec![unknown_file1, unknown_file2]);
            // Demultiplexing
            let (stats, unk_status) = demux::pe_demux(
                forward,
                reverse,
                format,
                utils::to_niffler_level(raw_level),
                &barcode_info,
                mismatch,
                &mut nb_records,
            )?;
            if !quiet {
                for (key, value) in stats.iter() {
                    info!(
                        "{} records found for {} barcode",
                        value,
                        String::from_utf8_lossy(key)
                    );
                }
            }
            if unk_status == *"truetrue" {
                fs::remove_file(future_unk_path1)?;
                fs::remove_file(future_unk_path2)?;
            } else if unk_status == *"falsetrue" {
                fs::remove_file(future_unk_path2)?;
            } else if unk_status == *"truefalse" {
                fs::remove_file(future_unk_path1)?;
            }
        }
    }

    if !quiet {
        // Finishing
        let duration = startime.elapsed();
        let seconds = duration.as_secs() % 60;
        let minutes = (duration.as_secs() / 60) % 60;
        let hours = (duration.as_secs() / 60) / 60;

        info!("Results are available in {}", output.display());
        info!("Walltime: {}h:{}m:{}s", hours, minutes, seconds,);
        info!("Used memory: {} bytes", sys.used_memory());
        info!("Thanks. Share. Come again!");
    }

    Ok(())
}
