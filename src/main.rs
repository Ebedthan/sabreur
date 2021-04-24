// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate bio;
extern crate clap;

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process;
use std::time::Instant;

use bio::io;
use clap::{App, Arg};

mod utils;

fn main() {
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

    utils::msg("sabreur v0.1.0 starting up!", quiet);
    if reverse.is_empty() {
        utils::msg("You are in single-end mode", quiet);
    } else {
        utils::msg("You are in paired-end mode", quiet);
    }

    // Handle output dir
    let output_exists = Path::new(output).exists();
    if output_exists && !force {
        utils::err(format!("Specified output folder: {}, already exists! Please change it using --out option or use --force to overwrite it.", output).as_str());
        process::exit(1);
    } else if output_exists && force {
        utils::msg(format!("Reusing directory {}", output).as_str(), quiet);
        fs::remove_dir_all(Path::new(output)).expect("Cannot remove existing directory");
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    } else if !output_exists {
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    }

    // Check File type: fasta or fastq
    let forward_file_ext = utils::get_file_type(forward);
    let reverse_file_ext = utils::get_file_type(reverse);

    if !reverse.is_empty() && forward_file_ext != reverse_file_ext {
        utils::err("Mismatched type of file supplied: one is fasta while other is fastq");
        process::exit(1);
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data = utils::read_file_to_string(barcode).unwrap();
    let barcode_fields = utils::split_line_by_tab(&barcode_data);

    if !utils::is_tab_delimited(&barcode_fields) {
        utils::err("Provided barcode file is not correcly formated");
        process::exit(1);
    }

    // Get forward file reader
    let (forward_reader, forward_compression) =
        utils::read_file(&Path::new(forward)).expect("Cannot read input file");

    if forward_compression == niffler::compression::Format::Gzip {
        utils::msg("Provided files are gzipped", quiet);
    }

    if mismatch != 0 {
        utils::warn(
            format!("You allowed {} mismatch in your barcode sequence", mismatch).as_str(),
            quiet,
        );
    }

    // Main processing of reads
    match forward_file_ext {
        Some(utils::FileType::Fasta) => match reverse.is_empty() {
            // single-end fasta mode
            true => {
                for b_vec in barcode_fields.iter() {
                    barcode_info.insert(b_vec[0], vec![b_vec[1]]);
                }
                let mut fa_forward_records = io::fasta::Reader::new(forward_reader).records();
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
                let mut fa_forward_records = io::fasta::Reader::new(forward_reader).records();
                let (reverse_reader, reverse_compression) =
                    utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                let mut fa_reverse_records = io::fasta::Reader::new(reverse_reader).records();
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
                let mut fq_forward_records = io::fastq::Reader::new(forward_reader).records();
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
                let mut fq_forward_records = io::fastq::Reader::new(forward_reader).records();
                let (reverse_reader, reverse_compression) =
                    utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                let mut fq_reverse_records = io::fastq::Reader::new(reverse_reader).records();

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
            utils::err("One of the provided file is not fasta nor fastq");
            process::exit(1);
        }
    }

    // Finishing
    let duration = startime.elapsed();
    let seconds = duration.as_secs() % 60;
    let minutes = (duration.as_secs() / 60) % 60;
    let hours = (duration.as_secs() / 60) / 60;

    utils::msg(
        format!("{} {}", "Results are available in", output).as_str(),
        quiet,
    );
    utils::msg(
        format!("Walltime: {}h:{}m:{}s", hours, minutes, seconds).as_str(),
        quiet,
    );
    utils::msg("Thanks. Share. Come again!", quiet);
}
