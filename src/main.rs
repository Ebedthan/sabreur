// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate anyhow;
extern crate clap;
extern crate exitcode;
extern crate human_panic;

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use anyhow::{anyhow, Context, Result};
use clap::{App, Arg};
use human_panic::setup_panic;

mod error;
mod utils;

const VERSION: &str = "0.3.0";

fn main() -> Result<()> {
    setup_panic!();

    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version(format!("v{}", VERSION).as_str())
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Fast, reliable and handy barcode demultiplexing tool for fastx files")
        .arg(
            Arg::with_name("BARCODE")
                .help("Input barcode file.")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("FORWARD")
                .help("Input forward fasta or fastq file. Can be gz, xz or bz2 compressed.")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::with_name("REVERSE")
                .help("Input reverse fasta or fastq file. Can be gz, xz or bz2 compressed.")
                .index(3),
        )
        .arg(
            Arg::with_name("mismatch")
                .help("Maximum number of mismatches allowed in a barcode")
                .short("m")
                .long("mismatch")
                .value_name("N")
                .default_value("0"),
        )
        .arg(
            Arg::with_name("output")
                .help("Output folder")
                .short("o")
                .long("out")
                .value_name("FOLDER")
                .default_value("sabreur_out"),
        )
        .arg(
            Arg::with_name("format")
                .help("Set output files compression format.")
                .long("format")
                .short("f")
                .takes_value(true)
                .possible_values(&["gz", "xz", "bz2"]),
        )
        .arg(
            Arg::with_name("level")
                .help("Set the compression level")
                .long("level")
                .short("l")
                .takes_value(true)
                .possible_values(&["1", "2", "3", "4", "5", "6", "7", "8", "9"])
                .default_value("1"),
        )
        .arg(
            Arg::with_name("force")
                .help("Force reuse of output directory")
                .long("force")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("quiet")
                .help("Decrease program verbosity")
                .short("q")
                .long("quiet")
                .takes_value(false),
        )
        .get_matches();

    // START ----------------------------------------------------------------
    let startime = Instant::now();
    let stdout = io::stdout();
    let mut ohandle = stdout.lock();
    let stderr = io::stderr();
    let mut ehandle = stderr.lock();

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
    let quiet = matches.is_present("quiet");

    if !quiet {
        writeln!(ohandle, "[INFO] sabreur v{} starting up!", VERSION)?;
        if reverse.is_empty() {
            writeln!(ohandle, "[INFO] You are in single-end mode")?;
        } else {
            writeln!(ohandle, "[INFO] You are in paired-end mode")?;
        }
    }

    // Change file compression format here for files extension
    if format != niffler::compression::Format::No {
        forward_file_compression = format;
        if !quiet {
            writeln!(
                ohandle,
                "[INFO] Output files will be {} compressed",
                utils::to_compression_ext(forward_file_compression)
            )?;
        }
    }

    // Exit if files does not have same types
    if !reverse.is_empty()
        && forward_file_type != (utils::get_file_type(reverse)?).0
    {
        writeln!(ehandle, "[ERROR] Mismatched type of file supplied: one is fasta while the other is fastq")?;
        process::exit(exitcode::DATAERR);
    }

    // Handle output dir
    let outdir_exists = Path::new(output).exists();
    if outdir_exists && !force {
        writeln!(ehandle, "[ERROR] Specified output folder '{}', already exists!\n[ERROR] Please change folder name using --out or use --force.", output)?;
        process::exit(exitcode::CANTCREAT);
    } else if outdir_exists && force {
        if !quiet {
            writeln!(ohandle, "[INFO] Reusing directory {}", output)?;
        }
        fs::remove_dir_all(Path::new(output)).with_context(|| anyhow!("Could not remove folder '{}'. Do you have permission to remove this folder?", output))?;
        fs::create_dir(Path::new(output)).with_context(|| anyhow!("Could not create folder '{}'. Do you have permission to create this folder?", output))?;
    } else if !outdir_exists {
        fs::create_dir(Path::new(output))?;
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data = fs::read_to_string(barcode)?;
    let barcode_fields = utils::split_by_tab(&barcode_data).unwrap();

    if mismatch != 0 && !quiet {
        writeln!(ohandle, "[WARN] Barcode mismatch allowed: {}", mismatch)?;
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
                    let mut file_path = PathBuf::from("");
                    file_path.push(output);
                    file_path.push(format!("{}{}", b_vec[1], ext));

                    let file = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path)?;

                    barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
                }

                // Create unknown file
                let mut unk_path = PathBuf::from("");
                unk_path.push(output);
                unk_path.push(format!("{}{}", "unknown.fa", ext));
                let future_unk_path = unk_path.clone();

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)?;

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                writeln!(ohandle, "[INFO] Demultiplexing ...")?;
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
                        writeln!(
                            ohandle,
                            "[INFO] {} records found for {} barcode",
                            value,
                            String::from_utf8_lossy(key)
                        )?;
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
                    let mut file_path = PathBuf::from("");
                    file_path.push(output);
                    file_path.push(format!("{}{}", b_vec[1], f_ext));

                    let mut file_path1 = PathBuf::from("");
                    file_path1.push(output);
                    file_path1.push(format!("{}{}", b_vec[2], r_ext));

                    let file1 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path)?;

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)?;

                    barcode_info
                        .insert(b_vec[0].as_bytes(), vec![file1, file2]);
                }

                // Create unknown files
                let mut unk_path1 = PathBuf::from("");
                unk_path1.push(output);
                unk_path1.push(format!("{}{}", "unknown_R1.fa", f_ext));

                let mut unk_path2 = PathBuf::from("");
                unk_path2.push(output);
                unk_path2.push(format!("{}{}", "unknown_R2.fa", r_ext));

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
                writeln!(ohandle, "[INFO] Demultiplexing ...")?;
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
                        writeln!(
                            ohandle,
                            "[INFO] {} records found for {} barcode",
                            value,
                            String::from_utf8_lossy(key)
                        )?;
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
                    let mut file_path = PathBuf::from("");
                    file_path.push(output);
                    file_path.push(format!("{}{}", b_vec[1], ext));

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
                unk_path.push(format!("{}{}", "unknown.fa", ext));

                let future_unk_path = unk_path.clone();

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)?;

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                writeln!(ohandle, "[INFO] Demultiplexing ...")?;
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
                        writeln!(
                            ohandle,
                            "[INFO] {} records found for {} barcode",
                            value,
                            String::from_utf8_lossy(key)
                        )?;
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
                    let mut file_path = PathBuf::from("");
                    file_path.push(output);
                    file_path.push(format!("{}{}", b_vec[1], f_ext));

                    let mut file_path1 = PathBuf::from("");
                    file_path1.push(output);
                    file_path1.push(format!("{}{}", b_vec[2], r_ext));

                    let file1 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path)?;

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)?;

                    barcode_info
                        .insert(b_vec[0].as_bytes(), vec![file1, file2]);
                }

                // Create unknown files
                let mut unk_path1 = PathBuf::from("");
                unk_path1.push(output);
                unk_path1.push(format!("{}{}", "unknown_R1.fq", f_ext));

                let mut unk_path2 = PathBuf::from("");
                unk_path2.push(output);
                unk_path2.push(format!("{}{}", "unknown_R2.fq", r_ext));

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
                writeln!(ohandle, "[INFO] Demultiplexing ...")?;
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
                        writeln!(
                            ohandle,
                            "[INFO] {} records found for {} barcode",
                            value,
                            String::from_utf8_lossy(key)
                        )?;
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
            writeln!(
                ehandle,
                "[ERROR] Supplied files are not fasta or fastq. Is there fas, fa, fasta, fq or fastq in the filename?"
            )?;
            process::exit(exitcode::DATAERR);
        }
    }

    if !quiet {
        // Finishing
        let duration = startime.elapsed();
        let seconds = duration.as_secs() % 60;
        let minutes = (duration.as_secs() / 60) % 60;
        let hours = (duration.as_secs() / 60) / 60;

        writeln!(ohandle, "[INFO] Results are available in {}", output)?;
        writeln!(
            ohandle,
            "[INFO] Walltime: {}h:{}m:{}s",
            hours, minutes, seconds,
        )?;
        writeln!(ohandle, "[INFO] Thanks. Share. Come again!")?;
    }

    Ok(())
}
