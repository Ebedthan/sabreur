// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate bio;
extern crate clap;
extern crate exitcode;
extern crate human_panic;

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process;
use std::time::Instant;

use bio::io::{fasta, fastq};
use clap::{App, Arg};
use human_panic::setup_panic;

mod utils;

const VERSION: &str = "0.1.2";

fn main() {
    setup_panic!();

    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version("v0.1.0")
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
    let forward = matches.value_of("FORWARD").unwrap();

    let mut reverse = "";
    if matches.is_present("REVERSE") {
        reverse = matches.value_of("REVERSE").unwrap();
    }

    let barcode = matches.value_of("BARCODE").unwrap();
    let output = matches.value_of("output").unwrap();
    let mis = matches.value_of("mismatch").unwrap().to_string();
    let mismatch = mis.parse::<i32>().unwrap();
    let force = matches.is_present("force");
    let quiet = matches.is_present("quiet");

    if !quiet {
        writeln!(io::stdout(), "[INFO] sabreur v{} starting up!", VERSION)
            .expect("Cannot write to stdout");
        if reverse.is_empty() {
            writeln!(io::stdout(), "[INFO] You are in single-end mode")
                .expect("Cannot write to stdout");
        } else {
            writeln!(io::stdout(), "[INFO] You are in paired-end mode")
                .expect("Cannot write to stdout");
        }
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
        fs::remove_dir_all(Path::new(output))
            .expect("Cannot remove existing directory");
        fs::create_dir(Path::new(output))
            .expect("Cannot create output directory");
    } else if !output_exists {
        fs::create_dir(Path::new(output))
            .expect("Cannot create output directory");
    }

    // Check File type: fasta or fastq
    let forward_file_ext = utils::get_file_type(forward);
    let reverse_file_ext = utils::get_file_type(reverse);

    if !reverse.is_empty() && forward_file_ext != reverse_file_ext {
        writeln!(
            io::stderr(),
            "[ERROR] Mismatched type of file supplied: one is fasta while other is fastq"
        )
        .expect("Cannot write to stderr");
        process::exit(exitcode::DATAERR);
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

    // Main processing of reads
    match forward_file_ext {
        Some(utils::FileType::Fasta) => match reverse.is_empty() {
            // single-end fasta mode
            true => {
                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0].as_bytes(), vec![b_vec[1]]);
                }

                // Get forward reader and compression format
                let (forward_reader, compression) = utils::read_file(forward)
                    .expect("Cannot read input file");
                let mut fa_forward_records =
                    fasta::Reader::new(forward_reader).records();

                // Demultiplexing
                let stats = utils::se_fa_demux(
                    &mut fa_forward_records,
                    compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demutiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {} contains {} records",
                            key, value
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
            // paired-end fasta mode
            false => {
                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(
                        b_vec[0].as_bytes(),
                        vec![b_vec[1], b_vec[2]],
                    );
                }

                // Get files reader and compressions format
                let (forward_reader, compression) = utils::read_file(forward)
                    .expect("Cannot read input file");
                let mut fa_forward_records =
                    fasta::Reader::new(forward_reader).records();

                let (reverse_reader, _compression) =
                    utils::read_file(reverse)
                        .expect("Cannot read input file");
                let mut fa_reverse_records =
                    fasta::Reader::new(reverse_reader).records();

                // Demultiplexing
                let stats = utils::pe_fa_demux(
                    &mut fa_forward_records,
                    &mut fa_reverse_records,
                    compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {} contains {} records",
                            key, value
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
        },
        Some(utils::FileType::Fastq) => match reverse.is_empty() {
            // single-end fastq mode
            true => {
                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0].as_bytes(), vec![b_vec[1]]);
                }

                // Get files reader and compression format
                let (forward_reader, compression) = utils::read_file(forward)
                    .expect("Cannot read input file");
                let mut fq_forward_records =
                    fastq::Reader::new(forward_reader).records();

                // Demultiplexing
                let stats = utils::se_fq_demux(
                    &mut fq_forward_records,
                    compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {} contains {} records",
                            key, value
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
            // paired-end fastq mode
            false => {
                // Read barcode data
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(
                        b_vec[0].as_bytes(),
                        vec![b_vec[1], b_vec[2]],
                    );
                }

                // Get files readers and compression
                let (forward_reader, compression) = utils::read_file(forward)
                    .expect("Cannot read input file");
                let mut fq_forward_records =
                    fastq::Reader::new(forward_reader).records();

                let (reverse_reader, _compression) =
                    utils::read_file(reverse)
                        .expect("Cannot read input file");
                let mut fq_reverse_records =
                    fastq::Reader::new(reverse_reader).records();
                // Demultiplexing
                let stats = utils::pe_fq_demux(
                    &mut fq_forward_records,
                    &mut fq_reverse_records,
                    compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");

                if !quiet {
                    for (key, value) in stats.iter() {
                        writeln!(
                            io::stdout(),
                            "[INFO] {} contains {} records",
                            key, value
                        )
                        .expect("Cannot write to stdout");
                    }
                }
            }
        },
        None => {
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
            hours, minutes, seconds
        )
        .expect("Cannot write to stdout");
        writeln!(io::stdout(), "Thanks. Share. Come again!")
            .expect("Cannot write to stdout");
    }

    process::exit(exitcode::OK);
}
