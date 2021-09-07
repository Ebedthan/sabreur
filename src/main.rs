// Copyright 2021 Anicet Ebou.
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
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use anyhow::{anyhow, Context, Result};
use clap::crate_version;
use log::{error, info, warn};
use sysinfo::{System, SystemExt};

mod app;
mod error;
mod utils;

// TODO: Check if supplied barcode file for se or pe is properly
// formated before giving it to the demultiplexing function
fn main() -> Result<()> {
    let mut sys = System::new_all();
    sys.refresh_all();
    let startime = Instant::now();

    // Define command-line arguments ----------------------------------------
    let matches = app::build_app().get_matches_from(env::args_os());

    // START ----------------------------------------------------------------
    let stderr = io::stderr();
    let mut ehandle = stderr.lock();

    // is --quiet option specified by the user?
    let quiet = matches.is_present("quiet");
    utils::setup_logging(quiet)?; // Settting up logging

    // Read command-line arguments
    let forward = matches
        .value_of("FORWARD")
        .with_context(|| anyhow!("Could not find input forward file"))?;
    let (forward_file_type, mut forward_file_compression) =
        utils::get_file_type(forward).with_context(|| {
            error::Error::UnableToDetectFileFormat {
                filename: forward.to_string(),
            }
        })?;

    let mut reverse = "";
    if matches.is_present("REVERSE") {
        reverse = matches
            .value_of("REVERSE")
            .with_context(|| anyhow!("Could not find input reverse file"))?;
    }

    let barcode = matches
        .value_of("BARCODE")
        .with_context(|| anyhow!("Could not find barcode file"))?;

    if !Path::new(barcode).exists() {
        writeln!(ehandle, "[ERROR] Barcode file not found. Is the path correct? with correct set of permissions?")?;
        process::exit(exitcode::DATAERR);
    }
    let output = matches.value_of("output").unwrap();
    let mis = matches.value_of("mismatch").unwrap().to_string();
    let mismatch = mis.parse::<i32>()?;

    // If user force output to be compressed even if input is not
    // add option to change compression of output
    let mut format = niffler::compression::Format::No;
    if matches.is_present("format") {
        format = utils::to_niffler_format(matches.value_of("format").unwrap())
            .with_context(|| {
                anyhow!(
                    "Could not convert compression format to niffler format"
                )
            })?;
    }

    let lv = matches.value_of("level").unwrap().to_string();
    let raw_level = lv.parse::<i32>()?;
    let force = matches.is_present("force");

    // Exit if files does not have same types
    if !reverse.is_empty()
        && forward_file_type != (utils::get_file_type(reverse)?).0
    {
        writeln!(ehandle, "[ERROR] Mismatched type of file supplied: one is fasta while the other is fastq")?;
        process::exit(exitcode::DATAERR);
    }

    info!("sabreur v{} starting up!", crate_version!());
    if reverse.is_empty() {
        info!("You are in single-end mode");
    } else {
        info!("You are in paired-end mode");
    }

    // Change file compression format here for files extension
    if format != niffler::compression::Format::No {
        forward_file_compression = format;
        info!(
            "Output files will be {} compressed",
            utils::to_compression_ext(forward_file_compression)
        );
    }

    // Handle output dir
    let outdir_exists = Path::new(output).exists();
    if outdir_exists && !force {
        error!("Specified output folder '{}', already exists!\nPlease change folder name using --out or use --force.", output);
        process::exit(exitcode::CANTCREAT);
    } else if outdir_exists && force {
        info!("Reusing directory {}", output);
        fs::remove_dir_all(Path::new(output)).with_context(|| anyhow!("Could not remove folder '{}'. Do you have permission to remove this folder?", output))?;
        fs::create_dir(Path::new(output)).with_context(|| anyhow!("Could not create folder '{}'. Do you have permission to create this folder?", output))?;
    } else if !outdir_exists {
        fs::create_dir(Path::new(output))?;
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data = fs::read_to_string(barcode)?;
    let barcode_fields = utils::split_by_tab(&barcode_data).unwrap();

    if mismatch != 0 {
        warn!("Barcode mismatch allowed: {}", mismatch);
    }

    let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

    // Main processing of reads
    match forward_file_type {
        utils::FileType::Fasta => match reverse.is_empty() {
            // single-end fasta mode
            true => {
                let ext = utils::to_compression_ext(forward_file_compression);

                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    let file_path = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[1], ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path)?;

                    barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
                }

                // Create unknown file
                let unk_path = utils::create_relpath_from(
                    [output, format!("{}{}", "unknown.fa", ext).as_str()]
                        .to_vec(),
                )
                .unwrap();
                let future_unk_path = unk_path.clone();

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)?;

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                let (stats, is_unk_empty) = utils::se_fa_demux(
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
                let (_reverse_file_type, mut reverse_file_compression) =
                    utils::get_file_type(reverse).unwrap();

                if format != niffler::compression::Format::No {
                    reverse_file_compression = format;
                }

                let f_ext = utils::to_compression_ext(forward_file_compression);
                let r_ext = utils::to_compression_ext(reverse_file_compression);

                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    let file_path1 = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[1], f_ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file_path2 = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[2], r_ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file1 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)?;

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path2)?;

                    barcode_info
                        .insert(b_vec[0].as_bytes(), vec![file1, file2]);
                }

                // Create unknown files
                let unk_path1 = utils::create_relpath_from(
                    [output, format!("{}{}", "unknown_R1.fa", f_ext).as_str()]
                        .to_vec(),
                )
                .unwrap();
                let unk_path2 = utils::create_relpath_from(
                    [output, format!("{}{}", "unknown_R2.fa", r_ext).as_str()]
                        .to_vec(),
                )
                .unwrap();

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
                let (stats, unk_status) = utils::pe_fa_demux(
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
        },
        utils::FileType::Fastq => match reverse.is_empty() {
            // single-end fastq mode
            true => {
                let ext = utils::to_compression_ext(forward_file_compression);

                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    let file_path = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[1], ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path)
                        .expect("Cannot open file");

                    barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
                }

                // Create unknown file
                let mut unk_path = PathBuf::from("");
                unk_path.push(output);
                unk_path.push(format!("{}{}", "unknown.fq", ext));

                let future_unk_path = unk_path.clone();

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)?;

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                let (stats, is_unk_empty) = utils::se_fq_demux(
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
            // paired-end fastq mode
            false => {
                let (_reverse_file_type, mut reverse_file_compression) =
                    utils::get_file_type(reverse).unwrap();

                if format != niffler::compression::Format::No {
                    reverse_file_compression = format;
                }

                let f_ext = utils::to_compression_ext(forward_file_compression);
                let r_ext = utils::to_compression_ext(reverse_file_compression);

                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    let file_path1 = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[1], f_ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file_path2 = utils::create_relpath_from(
                        [output, format!("{}{}", b_vec[2], r_ext).as_str()]
                            .to_vec(),
                    )
                    .unwrap();

                    let file1 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)?;

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path2)?;

                    barcode_info
                        .insert(b_vec[0].as_bytes(), vec![file1, file2]);
                }

                // Create unknown files
                let unk_path1 = utils::create_relpath_from(
                    [output, format!("{}{}", "unknown_R1.fq", f_ext).as_str()]
                        .to_vec(),
                )
                .unwrap();
                let unk_path2 = utils::create_relpath_from(
                    [output, format!("{}{}", "unknown_R2.fq", r_ext).as_str()]
                        .to_vec(),
                )
                .unwrap();

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
                let (stats, unk_status) = utils::pe_fq_demux(
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
        },
        utils::FileType::None => {
            error!(
                "Supplied files are not fasta or fastq. Is there fas, fa, fasta, fq or fastq in the filename?"
            );
            process::exit(exitcode::DATAERR);
        }
    }

    if !quiet {
        // Finishing
        let duration = startime.elapsed();
        let seconds = duration.as_secs() % 60;
        let minutes = (duration.as_secs() / 60) % 60;
        let hours = (duration.as_secs() / 60) / 60;

        info!("Results are available in {}", output);
        info!("Walltime: {}h:{}m:{}s", hours, minutes, seconds,);
        info!("Used memory: {} KB", sys.used_memory());
        info!("Thanks. Share. Come again!");
    }

    Ok(())
}
