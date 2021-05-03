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

const VERSION: &str = "0.1.1";

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

    let stdout = io::stdout();
    let stderr = io::stderr();
    let mut out_handle = stdout.lock();
    let mut err_handle = stderr.lock();

    if !quiet {
        writeln!(out_handle, "[INFO] sabreur v{} starting up!", VERSION)
            .expect("Cannot write to stdout");
        if reverse.is_empty() {
            writeln!(out_handle, "[INFO] You are in single-end mode")
                .expect("Cannot write to stdout");
        } else {
            writeln!(out_handle, "[INFO] You are in paired-end mode")
                .expect("Cannot write to stdout");
        }
    }

    // Handle output dir
    let output_exists = Path::new(output).exists();
    if output_exists && !force {
        writeln!(err_handle, "[ERROR] Specified output folder: {}, already exists! Please change it using --out option or use --force to overwrite it.", output).expect("Cannot write to stderr");
        process::exit(exitcode::CANTCREAT);
    } else if output_exists && force {
        if !quiet {
            writeln!(out_handle, "[INFO] Reusing directory {}", output)
                .expect("Cannot write to stdout");
        }
        fs::remove_dir_all(Path::new(output)).expect("Cannot remove existing directory");
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    } else if !output_exists {
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    }

    // Check File type: fasta or fastq
    let forward_file_ext = utils::get_file_type(forward);
    let reverse_file_ext = utils::get_file_type(reverse);

    if !reverse.is_empty() && forward_file_ext != reverse_file_ext {
        writeln!(
            err_handle,
            "[ERROR] Mismatched type of file supplied: one is fasta while other is fastq"
        )
        .expect("Cannot write to stderr");
        process::exit(exitcode::DATAERR);
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data = utils::read_file_to_string(barcode).unwrap();
    let barcode_fields = utils::split_line_by_tab(&barcode_data);

    if !utils::is_tab_delimited(&barcode_fields) {
        writeln!(
            err_handle,
            "[ERROR] Provided barcode file is not correcly formated"
        )
        .expect("Cannot write to stderr");
        process::exit(exitcode::DATAERR);
    }

    // Get forward file reader
    let (forward_reader, forward_compression) =
        utils::read_file(&Path::new(forward)).expect("Cannot read input file");

    if forward_compression == niffler::compression::Format::Gzip && !quiet {
        writeln!(out_handle, "[INFO] Provided files are gzipped").expect("Cannot write to stdout");
    }

    if mismatch != 0 && !quiet {
        writeln!(
            out_handle,
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
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0], vec![b_vec[1]]);
                }
                let mut fa_forward_records = fasta::Reader::new(forward_reader).records();
                utils::se_fa_demux(
                    &mut fa_forward_records,
                    forward_compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demutiplex file");
            }
            // paired-end fasta mode
            false => {
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0], vec![b_vec[1], b_vec[2]]);
                }
                let mut fa_forward_records = fasta::Reader::new(forward_reader).records();
                let (reverse_reader, reverse_compression) =
                    utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                let mut fa_reverse_records = fasta::Reader::new(reverse_reader).records();
                utils::pe_fa_demux(
                    &mut fa_forward_records,
                    &mut fa_reverse_records,
                    reverse_compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");
            }
        },
        Some(utils::FileType::Fastq) => match reverse.is_empty() {
            // single-end fastq mode
            true => {
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0], vec![b_vec[1]]);
                }
                let mut fq_forward_records = fastq::Reader::new(forward_reader).records();
                utils::se_fq_demux(
                    &mut fq_forward_records,
                    forward_compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");
            }
            // paired-end fastq mode
            false => {
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0], vec![b_vec[1], b_vec[2]]);
                }
                let mut fq_forward_records = fastq::Reader::new(forward_reader).records();
                let (reverse_reader, reverse_compression) =
                    utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                let mut fq_reverse_records = fastq::Reader::new(reverse_reader).records();

                utils::pe_fq_demux(
                    &mut fq_forward_records,
                    &mut fq_reverse_records,
                    reverse_compression,
                    &barcode_info,
                    mismatch,
                    output,
                )
                .expect("Cannot demultiplex file");
            }
        },
        None => {
            writeln!(
                err_handle,
                "[ERROR] One of the provided file is not fasta nor fastq"
            )
            .expect("Cannot write to stderr");
            process::exit(exitcode::DATAERR);
        }
    }

    // Finishing
    let duration = startime.elapsed();
    let seconds = duration.as_secs() % 60;
    let minutes = (duration.as_secs() / 60) % 60;
    let hours = (duration.as_secs() / 60) / 60;

    writeln!(out_handle, "[INFO] Results are available in {}", output)
        .expect("Cannot write to stdout");
    writeln!(
        out_handle,
        "[INFO] Walltime: {}h:{}m:{}s",
        hours, minutes, seconds
    )
    .expect("Cannot write to stdout");
    writeln!(out_handle, "Thanks. Share. Come again!").expect("Cannot write to stdout");

    process::exit(exitcode::OK);
}
