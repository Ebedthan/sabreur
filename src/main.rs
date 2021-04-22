// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate bio;
extern crate clap;
extern crate num_cpus;

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process;

use bio::io;
use clap::{App, Arg};
use log::{info, error};

mod utils;

fn main() {
    
    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version("v0.1.0")
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("A barcode demultiplexing tool")
        .arg(Arg::with_name("BARCODE")
            .help("Input barcode file [required]")
            .required(true)
            .index(1))
        .arg(Arg::with_name("FORWARD")
            .help("Input forward fasta or fastq file. Can be gzipped [required]")
            .validator(utils::is_fastx)
            .required(true)
            .index(2))
        .arg(Arg::with_name("REVERSE")
            .help("Input reverse fasta or fastq file. Can be gzipped")
            .validator(utils::is_fastx)
            .index(3))
        .arg(Arg::with_name("output")
            .help("Output folder")
            .short("o")
            .long("out")
            .value_name("FOLDER")
            .default_value("sabreur_out"))
        .arg(Arg::with_name("cpus")
            .help("Specify the number of threads")
            .short("c")
            .long("cpus")
            .default_value("1"))
        .arg(Arg::with_name("force")
            .help("Force reuse of output directory")
            .long("force")
            .takes_value(false))
        .arg(Arg::with_name("quiet")
            .help("Decrease program verbosity")
            .short("q")
            .long("quiet")
            .takes_value(false))
        .get_matches();

    // START ----------------------------------------------------------------
    // Handle verbosity setting
    let quiet = matches.is_present("quiet");
    utils::setup_logging(quiet).expect("Failed to initialize logging.");

    // Handle cpus 
    let acpus = matches.value_of("cpus").unwrap();
    let mut cpus = acpus.parse::<u8>().unwrap();
    let n_cpus = num_cpus::get();

    if cpus == 0 {
        cpus = n_cpus as u8;
    } else if cpus > n_cpus as u8 {
        cpus = n_cpus as u8;
    } else {
        cpus = cpus;
    }

    // Read args
    let forward = matches.value_of("FORWARD").unwrap();
    println!("{}", cpus);

    let mut reverse = "";
    if matches.is_present("REVERSE") {
        reverse = matches.value_of("REVERSE").unwrap();
    }

    let barcode = matches.value_of("BARCODE").unwrap();
    let output = matches.value_of("output").unwrap();
    let force = matches.is_present("force");

    // Handle output dir
    if Path::new(output).exists() && !force {
        error!("Specified output folder already exists! Please change it using --out option or use --force to overwrite it.");
        process::exit(1);
    } else if Path::new(output).exists() && force {
        info!("Reusing directory");
        fs::remove_dir(Path::new(output)).expect("Cannot remove existing directory");
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    } else if !Path::new(output).exists() {
        fs::create_dir(Path::new(output)).expect("Cannot create output directory");
    }

    // Check File type: fasta or fastq
    let forward_file_ext = utils::get_file_type(forward).unwrap();
    let reverse_file_ext = utils::get_file_type(reverse).unwrap();

    if reverse != "" && forward_file_ext != reverse_file_ext {
        error!("Mismatched type of file supplied: one is fasta while other is fastq");
        process::exit(1);
    }

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_data = utils::read_file_to_string(barcode).unwrap();
    let barcode_fields = utils::split_line_by_tab(&barcode_data);

    for b_vec in barcode_fields.iter() {
        if b_vec.len() == 2 {
            barcode_info.insert(
                b_vec[0],
                vec![b_vec[1]]
            );
        } else if b_vec.len() == 3 {
            barcode_info.insert(
                b_vec[0],
                vec![b_vec[1], b_vec[2]]
            );
        }

    }

    // Get files reader
    let (forward_reader, _forward_compression) = utils::read_file(&Path::new(forward)).expect("Cannot read input file");

    match forward_file_ext {
        utils::FileType::Fasta =>
            match reverse.is_empty() {
                // single-end fasta mode
                true => {
                    let mut fa_forward_records = io::fasta::Reader::new(forward_reader).records();
                    utils::se_fa_demux(&mut fa_forward_records,
                           &barcode_info,
                           output)
                          .expect("Cannot demutiplex file");
                },
                // paired-end fasta mode
                false => {
                    let mut fa_forward_records = io::fasta::Reader::new(forward_reader).records();
                    let (reverse_reader, _reverse_compression) = utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                    let mut fa_reverse_records = io::fasta::Reader::new(reverse_reader).records();
                    utils::pe_fa_demux(&mut fa_forward_records,
                                       &mut fa_reverse_records,
                                       &barcode_info,
                                       output)
                                      .expect("Cannot demultiplex file");
                }
            },
        utils::FileType::Fastq =>
            match reverse.is_empty() {
                // single-end fastq mode
                true => {
                    let mut fq_forward_records = io::fastq::Reader::new(forward_reader).records();
                    utils::se_fq_demux(&mut fq_forward_records,
                                       &barcode_info,
                                       output)
                                      .expect("Cannot demultiplex file");
                },
                // paired-end fastq mode
                false => {
                    let mut fq_forward_records = io::fastq::Reader::new(forward_reader).records();
                    let (reverse_reader, _reverse_compression) = utils::read_file(&Path::new(reverse)).expect("Cannot read input file");
                    let mut fq_reverse_records = io::fastq::Reader::new(reverse_reader).records();

                    utils::pe_fq_demux(&mut fq_forward_records,
                                       &mut fq_reverse_records, 
                                       &barcode_info, 
                                       output)
                                      .expect("Cannot demultiplex file");
                }
            },
    }

    info!("{}", format_args!("{} {}", "Done! Results are available in ", output));

}