// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::prelude::*;
use std::path::Path;
use std::process;

extern crate chrono;
extern crate fern;
extern crate log;
extern crate niffler;

use anyhow::Result;
use bio::io::{fasta, fastq};
use bio::pattern_matching::kmp::KMP;
use regex::Regex;

// setup_logging ------------------------------------------------------------

/// Set up the program logging using fern crate.
/// 
/// # Example
/// ```rust
/// let quiet: bool = true;
/// setup_logging(quiet).expect("Cannot initialize program logging")
/// ```
///
pub fn setup_logging(quiet: bool) -> Result<(), fern::InitError> {  
    let mut base_config = fern::Dispatch::new();

    base_config = match quiet {
        true => base_config.level(log::LevelFilter::Info)
                        .level_for("overly-verbose-target", log::LevelFilter::Warn),
        _ => base_config.level(log::LevelFilter::Trace),
    };

    let file_config = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.target(),
                record.level(),
                message
            ))
        })
        .chain(fern::log_file("maph.log")?);
    
        let stdout_config = fern::Dispatch::new()
        .format(|out, message, record| {
            // special format for debug messages coming from maph.
            if record.level() > log::LevelFilter::Info && record.target() == "maph" {
                out.finish(format_args!(
                    "---\nDEBUG: {}: {}\n---",
                    chrono::Local::now().format("%H:%M:%S"),
                    message
                ))
            } else {
                out.finish(format_args!(
                    "[{}][{}] {}",
                    chrono::Local::now().format("%H:%M:%S"),
                    record.level(),
                    message
                ))
            }
        })
        .chain(io::stdout());

    base_config
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

// read_file function -------------------------------------------------------

/// Get reader and compression format of file
/// 
/// # Example
/// ```rust
/// # use std::path::Path;
/// 
/// let path = Path::new("path/to/file");
/// let (reader, compression) = read_file(&path);
/// ```
/// 
pub fn read_file(path: &Path) -> Result<(Box<dyn io::Read>, niffler::compression::Format), > {
    let raw_in = Box::new(io::BufReader::new(
        File::open(path)?
    ));

    let (reader, compression) = niffler::get_reader(raw_in)
        .expect("Cannot read input fasta file");

    match compression {
        niffler::compression::Format::Gzip => Ok((reader, compression)),
        _ => { 
            eprintln!("Provided file is not gzipped.");
            process::exit(1)
        },
    }

}

// is_fastx function --------------------------------------------------------

/// Validate file type from file extension
/// 
/// # Example
/// ```rust
/// let filename = "myfile.fa.gz".to_string();
/// assert!(is_fastx(filename).is_ok())
/// ```
/// 
pub fn is_fastx(filename: String) -> Result<(), String> {
    let ext = vec!["fa", "fas", "fasta", "fastq", "fq", "gz"];
    
    let path = Path::new(&filename);
    let f_ext = path.extension().unwrap();
    let f_last = f_ext.to_str().unwrap();

    if ext.contains(&f_last){
        return Ok(());
    } else {
        return Err("Input file is not fasta nor fastq formatted".to_string());
    }
}

// get_file_type function ---------------------------------------------------

// FileType structure
#[derive(Debug, PartialEq)]
pub enum FileType {
    Fasta,
    Fastq,
}

/// Get file type from filename
/// 
/// # Example
/// ```rust
/// let filename = "myfile.fq";
/// let file_type = get_file_type(filename);
/// ```
/// 
pub fn get_file_type(filename: &str) -> Option<FileType> {
    if filename.contains(".fastq") || filename.contains(".fq") {
        Some(FileType::Fastq)
    } else if filename.contains(".fasta") || filename.contains(".fa") || filename.contains(".fas") {
        Some(FileType::Fasta)
    } else {
        None
    }
}


// read_file_to_string function ---------------------------------------------

/// Read a file into a string
/// 
/// # Example
/// ```rust
/// let filename = "myfile.txt";
/// let file_as_string = read_file_to_string(filename);
/// ```
/// 
pub fn read_file_to_string(filename: &str) -> io::Result<String> {
    let mut file = File::open(filename)?;
    let mut s = String::new();
    file.read_to_string(&mut s)?;
    Ok(s)
}

// split_line_by_tab function ------------------------------------------------------

/// Split a &str at each \t
/// 
/// # Example
/// ```rust
/// let mystring = "hello\tworld";
/// let string_fields = split_line_by_tab(mystring);
/// ```
/// 
pub fn split_line_by_tab<'a>(string: &'a str) -> Vec<Vec<&'a str>> {
    let tab_re = Regex::new(r"\t").unwrap();
    string.lines().map(|line| {
        tab_re.split(line).collect()
    }).collect()
}


// Barcode type -------------------------------------------------------------
pub type Barcode<'a> = HashMap<&'a str, Vec<&'a str>>;


// write_to_fa function -----------------------------------------------------

/// Write to provided data to a fasta file in append mode
/// 
/// # Example
/// ```rust
/// 
/// # use bio::io::fasta;
/// let filename = "myfile.fa";
/// let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
/// write_to_fa(filename, &record);
/// ```
/// 
pub fn write_to_fa<'a>(filename: &'a str, record: &'a fasta::Record) -> Result<()>{
    let file = OpenOptions::new()
                .append(true)
                .create(true)
                .open(filename)?;
    let handle = io::BufWriter::new(file);
    let mut writer = fasta::Writer::new(handle);
    let _write_res = writer.write_record(&record)
                        .expect("Cannot write to fasta file");

    Ok(())
}

/// Write to provided data to a fastq file in append mode
/// 
/// # Example
/// ```rust
/// 
/// # use bio::io::fastq;
/// let filename = "myfile.fq";
/// let record = fastq::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
/// write_to_fq(filename, &record);
/// ```
/// 
/// 
pub fn write_to_fq<'a>(filename: &'a str, record: &'a fastq::Record) -> Result<()>{
    let file = OpenOptions::new()
                .append(true)
                .create(true)
                .open(filename)?;
    let handle = io::BufWriter::new(file);
    let mut writer = fastq::Writer::new(handle);
    let _write_res = writer.write_record(&record)
                        .expect("Cannot write to fasta file");

    Ok(())
}

// se_fa_demux function --------------------------------------------------------
pub fn se_fa_demux(
    reader: fasta::Reader<std::boxed::Box<dyn std::io::Read>>, 
    barcode_data: &Barcode,
    unknown_file: &str) -> Result<()> {
        
        for record in reader.records() {
            let record = record.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(record.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fa(value[0], &record).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fa(unknown_file, &record).expect("Cannot write to unknown file");
            }
        }

        Ok(())
}

// se_fq_demux function --------------------------------------------------------
pub fn se_fq_demux(
    reader: fastq::Reader<std::boxed::Box<dyn std::io::Read>>, 
    barcode_data: &Barcode,
    unknown_file: &str) -> Result<()> {
        
        for record in reader.records() {
            let record = record.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(record.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fq(value[0], &record).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fq(unknown_file, &record).expect("Cannot write to unknown file");
            }
        }

        Ok(())
}

// pe_fa_demux function --------------------------------------------------------
pub fn pe_fa_demux(
    forward_reader: fasta::Reader<std::boxed::Box<dyn std::io::Read>>,
    reverse_reader: fasta::Reader<std::boxed::Box<dyn std::io::Read>>, 
    barcode_data: &Barcode,
    unknown_file: &str) -> Result<()> {
        
        for f_rec in forward_reader.records() {
            let f_rec = f_rec.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(f_rec.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fa(value[0], &f_rec).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fa(unknown_file, &f_rec).expect("Cannot write to unknown file");
            }
        }

        for r_rec in reverse_reader.records() {
            let r_rec = r_rec.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(r_rec.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fa(value[0], &r_rec).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fa(unknown_file, &r_rec).expect("Cannot write to unknown file");
            }
        }

        Ok(())
}

// pe_fq_demux function --------------------------------------------------------
pub fn pe_fq_demux(
    forward_reader: fastq::Reader<std::boxed::Box<dyn std::io::Read>>,
    reverse_reader: fastq::Reader<std::boxed::Box<dyn std::io::Read>>, 
    barcode_data: &Barcode,
    unknown_file: &str) -> Result<()> {
        
        for f_rec in forward_reader.records() {
            let f_rec = f_rec.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(f_rec.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fq(value[0], &f_rec).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fq(unknown_file, &f_rec).expect("Cannot write to unknown file");
            }
        }

        for r_rec in reverse_reader.records() {
            let r_rec = r_rec.unwrap();
            let mut matched = false;
            for (key, value) in barcode_data {
                let kmp = KMP::new(key.as_bytes());
                let occ: Vec<usize> = kmp.find_all(r_rec.seq().to_vec()).collect();
                if occ.contains(&0) {
                    write_to_fq(value[0], &r_rec).expect("Cannot write to output file");
                    matched = true;
                }
            }
            if matched == false {
                write_to_fq(unknown_file, &r_rec).expect("Cannot write to unknown file");
            }
        }

        Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_se_fa_demux() {
        let p = Path::new("tests/test2.fa.gz");
        let (fr, _cmp) = read_file(&p).expect("Cannot open");
        let reader = fasta::Reader::new(fr);
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("ACCGTA", vec!["tests/id1.fa"]);
        bc_data.insert("ATTGTT", vec!["tests/id2.fa"]);

        assert!(se_fa_demux(reader, &bc_data, "tests/unknown.fa").is_ok());

        std::fs::remove_file("tests/id1.fa").expect("Cannot delete tmp file");
        std::fs::remove_file("tests/id2.fa").expect("Cannot delete tmp file");
        std::fs::remove_file("tests/unknown.fa").expect("Cannot delete tmp file");
    }

    #[test]
    fn test_write_to_fa_is_ok() {
        let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");

        assert!((write_to_fa("tests/mytmp.fa", &record)).is_ok());

        let mut fa_records = fasta::Reader::from_file("tests/mytmp.fa")
                                .expect("Cannot read file.")
                                .records();

        while let Some(Ok(rec)) = fa_records.next() {
            assert_eq!(rec.id(), "id_str");
            assert_eq!(rec.desc(), Some("desc"));
            assert_eq!(rec.seq().to_vec(), b"ATCGCCG");
        }
        std::fs::remove_file("tests/mytmp.fa")
                    .expect("Cannot remove tmp file");
    }

    #[test]
    fn test_read_file() {
        let (mut reader, compression) = read_file(Path::new("tests/test.fa.gz")).unwrap();
        let mut contents = String::new();
        reader.read_to_string(&mut contents).expect("Error during file reading");

        assert_eq!(compression, niffler::compression::Format::Gzip);
        assert_eq!(contents, ">seqID1 desc\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n");
    }

    #[test]
    fn test_read_file_content_is_ok() {
        let (reader, _compression) = read_file(Path::new("tests/test.fa.gz")).unwrap();
        let fa_records = bio::io::fasta::Reader::new(reader).records();

        for record in fa_records {
            let record = record.unwrap();
            assert_eq!(record.id(), "seqID1");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq().to_vec(), b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        }
        
    }

    #[test]
    fn test_is_fastx() {
        let f1 = "myfile.fa";
        let f2 = "myfile.fq";
        let f3 = "myfile.fa.gz";
        let f4 = "myfile.txt";
        let f5 = "myfile.txt.xz";

        assert!(is_fastx(f1.to_string()).is_ok());
        assert!(is_fastx(f2.to_string()).is_ok());
        assert!(is_fastx(f3.to_string()).is_ok());
        assert!(is_fastx(f4.to_string()).is_err());
        assert!(is_fastx(f5.to_string()).is_err());

    }

    #[test]
    fn test_get_file_type(){
        let f1 = "myfile.fa";
        let f2 = "myfile.fq";
        let f3 = "myfile.fa.gz";
        let f4 = "myfile.txt";
        let f5 = "myfile.txt.xz";

        assert_eq!(get_file_type(f1), Some(FileType::Fasta));
        assert_eq!(get_file_type(f2), Some(FileType::Fastq));
        assert_eq!(get_file_type(f3), Some(FileType::Fasta));
        assert_ne!(get_file_type(f4), Some(FileType::Fastq));
        assert_eq!(get_file_type(f5), None);
    }

    #[test]
    fn test_read_file_to_string(){
        let mut file = OpenOptions::new()
                            .create(true)
                            .write(true)
                            .open("tests/mytmp.txt")
                            .unwrap();
        if let Err(e) = writeln!(file, "Hello\tWorld\tEarth\nBrian\twas\tthere\n") {
            eprintln!("Could not write to file: {}", e);
        }
        assert_eq!(
            read_file_to_string("tests/mytmp.txt").unwrap().trim(),
            "Hello\tWorld\tEarth\nBrian\twas\tthere");
        std::fs::remove_file("tests/mytmp.txt").expect("Cannot delete file");
    }

    #[test]
    fn test_split_line_by_tab(){
        let mystring = "Hello\tMars\tLinux\nBrian\tfork";
        let fields = split_line_by_tab(mystring);
        let myvec = vec![vec!["Hello","Mars","Linux"], vec!["Brian", "fork"]];
        assert_eq!(fields, myvec);
        
    }
}