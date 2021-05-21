// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate clap;
extern crate exitcode;
extern crate human_panic;

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use clap::{App, Arg};
use human_panic::setup_panic;

mod utils;

const VERSION: &str = "0.3.0";

fn main() {
    setup_panic!();

    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version(format!("v{}", VERSION).as_str())
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("A barcode demultiplexing tool for fasta and fastq files")
        .arg(
            Arg::with_name("BARCODE")
                .help("Input barcode file [required]")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("FORWARD")
                .help("Input forward fasta or fastq file. Can be gzipped [required]")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::with_name("REVERSE")
                .help("Input reverse fasta or fastq file. Can be gzipped")
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

    // Read args
    // Get filetype and compression format of input files
    let forward = matches.value_of("FORWARD").unwrap();
    let (forward_file_type, mut forward_file_compression) =
        utils::get_file_type(forward).unwrap();

    let mut reverse = "";

    if matches.is_present("REVERSE") {
        reverse = matches.value_of("REVERSE").unwrap();
    }
    let (reverse_file_type, mut reverse_file_compression) =
        utils::get_file_type(reverse).unwrap();

    let barcode = matches.value_of("BARCODE").unwrap();
    let output = matches.value_of("output").unwrap();
    let mis = matches.value_of("mismatch").unwrap().to_string();
    let mismatch = mis.parse::<i32>().unwrap();

    // If user force output to be compressed even if input is not
    // add option to change compression of output
    let mut format = niffler::compression::Format::No;
    if matches.is_present("format") {
        format = utils::to_niffler_format(matches.value_of("format").unwrap());
    }

    // Change file compression format here for files extension
    if format != niffler::compression::Format::No {
        forward_file_compression = format;
        reverse_file_compression = format;
    }

    let lv = matches.value_of("level").unwrap().to_string();
    let raw_level = lv.parse::<i32>().unwrap();
    let force = matches.is_present("force");
    let quiet = matches.is_present("quiet");

    if !quiet {
        writeln!(
            io::stdout(),
            "{}",
            format!("[INFO] sabreur v{} starting up!", VERSION).as_str()
        )
        .expect("Cannot write to stdout");
        if reverse.is_empty() {
            writeln!(io::stdout(), "[INFO] You are in single-end mode")
                .expect("Cannot write to stdout");
        } else {
            writeln!(io::stdout(), "[INFO] You are in paired-end mode")
                .expect("Cannot write to stdout");
        }
    }

    // Exit if files does not have same types
    if !reverse.is_empty() && forward_file_type != reverse_file_type {
        writeln!(io::stderr(), "[ERROR] Mismatched type of file supplied: one is fasta while the other is fastq").expect("Cannot write to stderr");
        process::exit(exitcode::DATAERR);
    }

    // Handle output dir
    let output_exists = Path::new(output).exists();
    if output_exists && !force {
        writeln!(io::stderr(), "[ERROR] Specified output folder: {}, already exists!\nPlease change it using --out option or use --force to overwrite it.", output).expect("Cannot write to stderr");
        process::exit(exitcode::CANTCREAT);
    } else if output_exists && force {
        if !quiet {
            writeln!(io::stdout(), "[INFO] Reusing directory {}", output)
                .expect("Cannot write to stdout");
        }
        fs::remove_dir_all(Path::new(output)).expect("Cannot remove directory");
        fs::create_dir(Path::new(output)).expect("Cannot create directory");
    } else if !output_exists {
        fs::create_dir(Path::new(output)).expect("Cannot create directory");
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data =
        fs::read_to_string(barcode).expect("Cannot read barcode file");
    let barcode_fields = utils::split_by_tab(&barcode_data).unwrap();

    if mismatch != 0 && !quiet {
        writeln!(
            io::stdout(),
            "[WARN] You allowed {} mismatch in your barcode sequence",
            mismatch
        )
        .expect("Cannot write to stdout");
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
                        .open(file_path)
                        .expect("Cannot open file");

                    barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
                }

                // Create unknown file
                let mut unk_path = PathBuf::from("");
                unk_path.push(output);
                unk_path.push(format!("{}{}", "unknown.fa", ext));

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)
                    .expect("Cannot create file");

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                let stats = utils::se_fa_demux(
                    forward,
                    format,
                    utils::to_niffler_level(raw_level),
                    &barcode_info,
                    mismatch,
                    &mut nb_records,
                )
                .expect("Cannot demutiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {}",
                            format!(
                                "{} records found for {} barcode",
                                value,
                                String::from_utf8_lossy(key)
                            )
                            .as_str(),
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
            // paired-end fasta mode
            false => {
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
                        .open(file_path)
                        .expect("Cannot open file");

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)
                        .expect("Cannot open file");

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

                let unknown_file1 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path1)
                    .expect("Cannot create file");

                let unknown_file2 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path2)
                    .expect("Cannot create file");

                barcode_info.insert(b"XXX", vec![unknown_file1, unknown_file2]);

                // Demultiplexing
                let stats = utils::pe_fa_demux(
                    forward,
                    reverse,
                    format,
                    utils::to_niffler_level(raw_level),
                    &barcode_info,
                    mismatch,
                    &mut nb_records,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {}",
                            format!(
                                "{} records found for {} barcode",
                                value,
                                String::from_utf8_lossy(key)
                            )
                            .as_str(),
                        )
                        .expect("Cannot write to stdout");
                    }
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
                        .open(b_vec[1])
                        .expect("Cannot open file");

                    barcode_info.insert(b_vec[0].as_bytes(), vec![file]);
                }

                // Create unknown file
                let mut unk_path = PathBuf::from("");
                unk_path.push(output);
                unk_path.push(format!("{}{}", "unknown.fa", ext));

                let unknown_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path)
                    .expect("Cannot create file");

                barcode_info.insert(b"XXX", vec![unknown_file]);

                // Demultiplexing
                let stats = utils::se_fq_demux(
                    forward,
                    format,
                    utils::to_niffler_level(raw_level),
                    &barcode_info,
                    mismatch,
                    &mut nb_records,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {}",
                            format!(
                                "{} records found for {} barcode",
                                value,
                                String::from_utf8_lossy(key)
                            )
                            .as_str(),
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
            // paired-end fastq mode
            false => {
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
                        .open(file_path)
                        .expect("Cannot open file");

                    let file2 = fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(file_path1)
                        .expect("Cannot open file");

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

                let unknown_file1 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path1)
                    .expect("Cannot create file");

                let unknown_file2 = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(unk_path2)
                    .expect("Cannot create file");

                barcode_info.insert(b"XXX", vec![unknown_file1, unknown_file2]);

                // Demultiplexing
                let stats = utils::pe_fq_demux(
                    forward,
                    reverse,
                    format,
                    utils::to_niffler_level(raw_level),
                    &barcode_info,
                    mismatch,
                    &mut nb_records,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {} records found for {} barcode",
                            value,
                            String::from_utf8_lossy(key)
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
        },
        utils::FileType::None => {
            writeln!(
                io::stderr(),
                "[ERROR] One of the provided file is not fasta nor fastq"
            )
            .expect("Cannot write to stderr");
            process::exit(exitcode::DATAERR);
        }
    }

    if !quiet {
        // Finishing
        let duration = startime.elapsed();
        let seconds = duration.as_secs() % 60;
        let minutes = (duration.as_secs() / 60) % 60;
        let hours = (duration.as_secs() / 60) / 60;

        writeln!(io::stdout(), "[INFO] Results are available in {}", output)
            .expect("Cannot write to stdout");
        writeln!(
            io::stdout(),
            "[INFO] Walltime: {}h:{}m:{}s",
            hours,
            minutes,
            seconds,
        )
        .expect("Cannot write to stdout");
        writeln!(io::stdout(), "[INFO] Thanks. Share. Come again!")
            .expect("Cannot write to stdout");
    }

    process::exit(exitcode::OK);
}
