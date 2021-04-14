/*
Copyright (c) 2021 Anicet Ebou <anicet.ebou@gmail.com>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate fern;
extern crate log;
extern crate niffler;
extern crate relative_path;
extern crate num_cpus;
extern crate regex;

use clap::{App, Arg, SubCommand};
use log::info;
use std::collections::HashMap;
use std::io;

/* mod declaration */
mod utils;

fn main() -> io::Result<()> {
    
    // Define command-line arguments ----------------------------------------
    let matches = App::new("sabreur")
        .version("v0.1.0")
        .author("Anicet Ebou <anicet.ebou@gmail.com>")
        .about("A barcode demultiplexing tool")
        .subcommand(SubCommand::with_name("se")
                    .about("single-end mode")
                    .arg(Arg::with_name("file")
                        .help("Input fasta or fastq file. Can be gzipped.")
                        .short("f")
                        .long("file")
                        .value_name("FILE")
                        .validator(utils::is_fastx)
                        .required(true))
                    .arg(Arg::with_name("barcode")
                        .help("Input barcode file.")
                        .short("b")
                        .long("barcode")
                        .value_name("FILE")
                        .required(true))
                    .arg(Arg::with_name("unknown")
                        .help("Output file for sequence with no barcodes")
                        .short("u")
                        .long("unknown")
                        .value_name("FILE")
                        .required(true)))
        .subcommand(SubCommand::with_name("pe")
                        .about("paired-end mode")
                        .arg(Arg::with_name("forward")
                            .help("Input forward fasta or fastq file. Can be gzipped.")
                            .short("f")
                            .long("forward")
                            .value_name("FILE")
                            .required(true))
                        .arg(Arg::with_name("reverse")
                            .help("Input reverse fasta or fastq file. Can be gzipped.")
                            .short("r")
                            .long("reverse")
                            .value_name("FILE")
                            .required(true))
                        .arg(Arg::with_name("barcode")
                            .help("Input barcode file.")
                            .short("b")
                            .long("barcode")
                            .value_name("FILE")
                            .required(true))
                        .arg(Arg::with_name("unknown")
                            .help("Output file for sequence with no barcodes")
                            .short("u")
                            .long("unknown")
                            .value_name("FILE")
                            .required(true)))
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
    
    // Verbosity
    let quiet = matches.is_present("quiet");
    utils::setup_logging(quiet).expect("Failed to initialize logging.");

    // START ----------------------------------------------------------------
    info!("sabreur v0.1 starting up!");
    let cpus = matches.value_of("cpus").unwrap();
    if let Some(matches) = matches.subcommand_matches("se") {
        let infile = matches.value_of("file").unwrap();
        let barcode = matches.value_of("barcode").unwrap();
        let unknown = matches.value_of("unknown").unwrap();

        let mut barcode_info = HashMap::new();

        let barcode_reader = utils::read_barcode(barcode).unwrap();

        let barcode_fields = utils::split_line(&barcode_reader);
        
        for b_vec in barcode_fields.iter() {
            barcode_info.insert(utils::BarcodeOut::new(b_vec[2], b_vec[1]), b_vec[0]);
        }

    } else if let Some(matches) = matches.subcommand_matches("pe") {
        let infile = matches.value_of("forward").unwrap();
        let reverse = matches.value_of("reverse").unwrap();
        let barcode = matches.value_of("barcode").unwrap();
        let unknown = matches.value_of("unknown").unwrap();
        println!("{}", format!("infile: {}, reverse: {}, barcode: {}, unknown: {}, cpus: {}",
                    infile, reverse, barcode, unknown, cpus));
    }
    
    info!("Finished");

    Ok(())
}