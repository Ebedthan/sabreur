// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate fern;
extern crate log;
extern crate niffler;
extern crate relative_path;
extern crate num_cpus;
extern crate regex;


use clap::{App, Arg};
use log::{info, error};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};


mod utils;

fn main() {
    
    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version("v0.1.0")
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("A barcode demultiplexing tool")
        .arg(Arg::with_name("FORWARD")
            .help("Input forward fasta or fastq file. Can be gzipped.")
            .validator(utils::is_fastx)
            .required(true)
            .index(1))
        .arg(Arg::with_name("REVERSE")
            .help("Input reverse fasta or fastq file. Can be gzipped.")
            .validator(utils::is_fastx)
            .index(2))
        .arg(Arg::with_name("barcode")
            .help("Input barcode file.")
            .short("b")
            .long("barcode")
            .value_name("FILE")
            .required(true))
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
        .arg(Arg::with_name("quiet")
            .help("Decrease program verbosity.")
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

    let barcode = matches.value_of("barcode").unwrap();
    let output = matches.value_of("output").unwrap();

    // Check File type: fasta or fastq
    if reverse != "" && utils::get_file_type(forward).unwrap() != utils::get_file_type(reverse).unwrap() {
        error!("Mismatched type of file supplied: one is fasta while other is fastq");
    }

    // Create output directory
    fs::create_dir(Path::new(output)).expect("");

    // Define filename for unknown reads
    let file_ext = utils::get_file_type(forward).unwrap();

    let mut u_str = "";
    if file_ext == utils::FileType::Fasta {
        u_str = "unknown.fa"
    } else if file_ext == utils::FileType::Fastq {
        u_str = "unknown.fq"
    }
    let u_path: PathBuf = [output, u_str].iter().collect();
    let _unknown = fs::File::create(u_path);

    // Read data from barcode file
    let mut barcode_info: utils::Barcode = HashMap::new();
    let barcode_reader = utils::read_file_to_string(barcode).unwrap();
    let barcode_fields = utils::split_line_by_tab(&barcode_reader);

    for b_vec in barcode_fields.iter() {
        if b_vec.len() == 2 {
            barcode_info.insert(
                b_vec[0],
                utils::BarcodeOut {forward: b_vec[1], reverse: ""}
            );
        } else if b_vec.len() == 3 {
            barcode_info.insert(
                b_vec[0],
                utils::BarcodeOut {forward: b_vec[1], reverse: b_vec[2]}
            );
        }

    }

    info!("Finished");

}