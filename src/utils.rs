// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::Write;
use std::io::{self};
use std::path::PathBuf;
use std::process;

extern crate colored;
extern crate niffler;

use anyhow::{anyhow, Result};
use bio::io::{fasta, fastq};
use colored::*;

// msg function

/// Print given text to stdout
///
/// # Example
/// ```rust
/// msg("Hello World!");
/// ```
///
pub fn msg(msg: &str) {
    writeln!(io::stdout(), "{}", format!("[INFO] {}", msg))
        .expect("Cannot write to stdout");
}

// err function

/// Print given text to stderr in red
///
/// # Example
/// ```rust
/// err("This is an error");
/// ```
///
pub fn err(error: &str) {
    writeln!(io::stderr(), "{}", format!("[ERROR] {}", error).red())
        .expect("Cannot write to stderr");
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
pub fn read_file(
    filename: &str,
) -> Result<(Box<dyn io::Read>, niffler::compression::Format)> {
    let raw_in = Box::new(io::BufReader::new(File::open(filename)?));

    let (reader, compression) =
        niffler::get_reader(raw_in).expect("Cannot read input fasta file");

    match compression {
        niffler::compression::Format::Gzip => Ok((reader, compression)),
        niffler::compression::Format::No => Ok((reader, compression)),
        _ => {
            writeln!(io::stderr(), "[ERROR] Only gzipped files are supported")
                .expect("Cannot write to stderr");
            process::exit(1);
        }
    }
}

/// get_file_type function --------------------------------------------------

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
    } else if filename.contains(".fasta")
        || filename.contains(".fa")
        || filename.contains(".fas")
    {
        Some(FileType::Fasta)
    } else {
        None
    }
}

// split_line_by_tab function -----------------------------------------------

/// Split a &str at each \t
///
/// # Example
/// ```rust
/// let mystring = "hello\tworld";
/// let string_fields = split_line_by_tab(mystring);
/// ```
///
pub fn split_by_tab(string: &str) -> Result<Vec<Vec<&str>>> {
    if string.contains('\t') {
        Ok(string
            .lines()
            .map(|line| line.split('\t').collect())
            .collect())
    } else {
        Err(anyhow!("The barcode file is not tab-delimited"))
    }
}

// Barcode type -------------------------------------------------------------
pub type Barcode<'a> = HashMap<&'a [u8], Vec<std::fs::File>>;

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
pub fn write_to_fa<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fasta::Record,
) -> Result<()> {
    let handle = niffler::get_writer(
        Box::new(file),
        compression,
        niffler::compression::Level::One,
    )?;

    let mut writer = fasta::Writer::new(handle);
    let _write_res = writer.write_record(&record)?;

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
pub fn write_to_fq<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fastq::Record,
) -> Result<()> {
    let handle = niffler::get_writer(
        Box::new(file),
        compression,
        niffler::compression::Level::One,
    )?;

    let mut writer = fastq::Writer::new(handle);
    let _write_res = writer.write_record(&record)?;

    Ok(())
}

// bc_cmp fn ----------------------------------------------------------------

/// Compare provided barcode with a sequence
///
/// # Example
/// ```rust
///
/// let bc = "ATCT";
/// let seq = "ATCTGGGCCAAATTT";
/// bc_cmp(bc, seq);
///
/// ```
///
pub fn bc_cmp(bc: &[u8], seq: &[u8], mismatch: i32) -> bool {
    // This wonderful line below compute the number of
    // character mismatch between two strings
    let nb_mismatch: i32 = bc
        .iter()
        .zip(seq.iter())
        .map(|(a, b)| (a != b) as i32)
        .sum();

    nb_mismatch <= mismatch
}

// se_fa_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of single-end file
///
/// # Example
/// ```rust
///
/// let p = Path::new("file.fa.gz");
/// let out = "outfolder";
/// let (fr, cmp) = read_file(&p).expect("Cannot open");
/// let mut records = fasta::Reader::new(fr).records();
/// let mut bc_data: Barcode = HashMap::new();
/// bc_data.insert("ACCGTA", vec!["id1.fa"]);
/// bc_data.insert("ATTGTT", vec!["id2.fa"]);
///
/// assert!(se_fa_demux(&mut records, cmp, &bc_data, out).is_ok());
/// ```
///
///
pub fn se_fa_demux<'a>(
    forward_records: &'a mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &'a Barcode,
    mismatch: i32,
    out: &'a str,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<&'a mut HashMap<&'a [u8], i32>> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = ".gz";
    } else {
        ext = "";
    }

    let mut unk_path = PathBuf::from("");
    unk_path.push(out);
    unk_path.push(format!("{}{}", "unknown.fa", ext));
    let unknown_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path)
        .expect("Cannot create file");

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fa(
                    &barcode_data.get(i).unwrap()[0],
                    compression,
                    &f_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown.fq")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fa(&unknown_file, compression, &f_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    Ok(nb_records)
}

// se_fq_demux function -----------------------------------------------------

/// Demultiplex a fastq::Record of single-end file
///
/// # Example
/// ```rust
///
/// let p = Path::new("file.fq.gz");
/// let out = "outfolder";
/// let (fr, cmp) = read_file(&p).expect("Cannot open");
/// let mut records = fasta::Reader::new(fr).records();
/// let mut bc_data: Barcode = HashMap::new();
/// bc_data.insert("ACCGTA", vec!["id1.fa"]);
/// bc_data.insert("ATTGTT", vec!["id2.fa"]);
///
/// assert!(se_fq_demux(&mut records, cmp, &bc_data, out).is_ok());
/// ```
///
///
pub fn se_fq_demux<'a>(
    forward_records: &'a mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &'a Barcode,
    mismatch: i32,
    out: &'a str,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<&'a mut HashMap<&'a [u8], i32>> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = ".gz";
    } else {
        ext = "";
    }

    let mut unk_path = PathBuf::from("");
    unk_path.push(out);
    unk_path.push(format!("{}{}", "unknown.fq", ext));
    let unknown_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path)
        .expect("Cannot create file");

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fq(
                    &barcode_data.get(i).unwrap()[0],
                    compression,
                    &f_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown.fq")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fq(&unknown_file, compression, &f_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    Ok(nb_records)
}

// pe_fa_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fa_demux<'a>(
    forward_records: &'a mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    reverse_records: &'a mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &'a Barcode,
    mismatch: i32,
    out: &'a str,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<&'a mut HashMap<&'a [u8], i32>> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = ".gz";
    } else {
        ext = "";
    }

    let mut unk_path1 = PathBuf::from("");
    unk_path1.push(out);
    unk_path1.push(format!("{}{}", "unknown_R1.fa", ext));
    let mut unk_path2 = PathBuf::from("");
    unk_path2.push(out);
    unk_path2.push(format!("{}{}", "unknown_R2.fa", ext));
    let unknown_file1 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path1)
        .expect("Cannot create file");
    let unknown_file2 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path2)
        .expect("Cannot create file");

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fa(
                    &barcode_data.get(i).unwrap()[0],
                    compression,
                    &f_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown_R1.fa")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fa(&unknown_file1, compression, &f_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &r_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fa(
                    &barcode_data.get(i).unwrap()[1],
                    compression,
                    &r_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown_R2.fa")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fa(&unknown_file2, compression, &r_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    Ok(nb_records)
}

// pe_fq_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fq_demux<'a>(
    forward_records: &'a mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    reverse_records: &'a mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &'a Barcode,
    mismatch: i32,
    out: &'a str,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<&'a mut HashMap<&'a [u8], i32>> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = ".gz";
    } else {
        ext = "";
    }

    let mut unk_path1 = PathBuf::from("");
    unk_path1.push(out);
    unk_path1.push(format!("{}{}", "unknown_R1.fq", ext));
    let mut unk_path2 = PathBuf::from("");
    unk_path2.push(out);
    unk_path2.push(format!("{}{}", "unknown_R2.fq", ext));
    let unknown_file1 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path1)
        .expect("Cannot create file");
    let unknown_file2 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(unk_path2)
        .expect("Cannot create file");

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fq(
                    &barcode_data.get(i).unwrap()[0],
                    compression,
                    &f_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown_R1.fq")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fq(&unknown_file1, compression, &f_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut iter = my_vec.iter();
        let res = iter.find(|&&x| bc_cmp(x, &r_rec.seq()[..bc_len], mismatch));

        match res {
            Some(i) => {
                nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
                write_to_fq(
                    &barcode_data.get(i).unwrap()[1],
                    compression,
                    &r_rec,
                )
                .expect("Cannot write to output file");
            }
            None => {
                nb_records
                    .entry(b"unknown_R2.fq")
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                write_to_fq(&unknown_file2, compression, &r_rec)
                    .expect("Cannot write to unknown file");
            }
        }
    }

    Ok(nb_records)
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::prelude::*;

    // read_file tests ------------------------------------------------------
    #[test]
    fn test_read_file() {
        let (mut reader, compression) = read_file("tests/test.fa.gz").unwrap();
        let mut contents = String::new();
        reader
            .read_to_string(&mut contents)
            .expect("Error during file reading");

        assert_eq!(compression, niffler::compression::Format::Gzip);
        assert_eq!(contents, ">seqID1 desc\nATCGATCGATCGATC\n");
    }

    #[test]
    fn test_read_file_content_is_ok() {
        let (reader, _compression) = read_file("tests/test.fa.gz").unwrap();
        let fa_records = bio::io::fasta::Reader::new(reader).records();

        for record in fa_records {
            let record = record.unwrap();
            assert_eq!(record.id(), "seqID1");
            assert_eq!(record.desc(), Some("desc"));
            assert_eq!(record.seq().to_vec(), b"ATCGATCGATCGATC");
        }
    }

    // bc_cmp tests ---------------------------------------------------------
    #[test]
    fn test_bc_cmp_ok() {
        let seq = b"ATCGATCGATCG";
        let bc = b"ATCG";

        assert!(bc_cmp(bc, seq, 0));
    }

    #[test]
    fn test_bc_cmp_not_ok() {
        let bc = b"TGCA";
        let seq = b"ATCGATCGATCG";

        assert_eq!(bc_cmp(bc, seq, 0), false);
    }

    #[test]
    fn test_bc_cmp_mismatch_ok() {
        let bc = b"AACG";
        let seq = b"ATCGATCGATCG";

        assert!(bc_cmp(bc, seq, 1));
    }

    #[test]
    fn test_bc_cmp_mismatch_not_ok() {
        let bc = b"AACG";
        let seq = b"ATCGATCGATCG";

        assert_eq!(bc_cmp(bc, seq, 0), false);
    }

    // se_fa_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fa_demux() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fa.gz").expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fa_demux(
            &mut records,
            cmp,
            &bc_data,
            0,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m1() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fa.gz").expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fa_demux(
            &mut records,
            cmp,
            &bc_data,
            1,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m2() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fa.gz").expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fa_demux(
            &mut records,
            cmp,
            &bc_data,
            2,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    // se_fq_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fq_demux() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fq.gz").expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fq_demux(
            &mut records,
            cmp,
            &bc_data,
            0,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m1() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fq.gz").expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fq_demux(
            &mut records,
            cmp,
            &bc_data,
            1,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m2() {
        let out = "tests";
        let (fr, cmp) =
            read_file("tests/test2.fq.gz").expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();
        let file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id1.fa")
            .expect("cannot open file");
        let file2 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/id2.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);

        assert!(se_fq_demux(
            &mut records,
            cmp,
            &bc_data,
            2,
            out,
            &mut nb_records
        )
        .is_ok());
    }

    // write_to_fa tests ----------------------------------------------------
    #[test]
    fn test_write_to_fa_is_ok() {
        let record =
            fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
        let cmp = niffler::compression::Format::Gzip;
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/mytmp.fa")
            .expect("cannot open file");
        assert!((write_to_fa(&file, cmp, &record)).is_ok());

        let mut fa_records = fasta::Reader::from_file("tests/mytmp.fa")
            .expect("Cannot read file.")
            .records();

        while let Some(Ok(rec)) = fa_records.next() {
            assert_eq!(rec.id(), "id_str");
            assert_eq!(rec.desc(), Some("desc"));
            assert_eq!(rec.seq(), b"ATCGCCG");
        }
    }

    // write_to_fq tests ----------------------------------------------------
    #[test]
    fn test_write_to_fq_is_ok() {
        let record = fastq::Record::with_attrs(
            "id_str",
            Some("desc"),
            b"ATCGCCG",
            b"QQQQQQQ",
        );
        let cmp = niffler::compression::Format::Gzip;
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/mytmp.fa")
            .expect("cannot open file");
        assert!((write_to_fq(&file, cmp, &record)).is_ok());

        let mut fa_records = fastq::Reader::from_file("tests/mytmp.fq")
            .expect("Cannot read file.")
            .records();

        while let Some(Ok(rec)) = fa_records.next() {
            assert_eq!(rec.id(), "id_str");
            assert_eq!(rec.desc(), Some("desc"));
            assert_eq!(rec.seq(), b"ATCGCCG");
            assert_eq!(rec.qual(), b"QQQQQQQ");
        }
    }

    // get_file_type tests -------------------------------------------------
    #[test]
    fn test_get_file_type() {
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

    // split_by_tab tests ---------------------------------------------------
    #[test]
    fn test_split_by_tab() {
        let mystring = "Hello\tWorld\tEarth\nBrian\twas\tthere";
        let fields = split_by_tab(mystring).unwrap();
        assert_eq!(
            fields,
            [["Hello", "World", "Earth"], ["Brian", "was", "there"]]
        );
    }

    #[test]
    #[should_panic]
    fn test_split_by_tab_not_ok() {
        let mystring = "HelloWorldEarth\nBrianwasthere";
        split_by_tab(mystring).unwrap();
    }
}
