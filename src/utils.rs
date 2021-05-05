// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process;

extern crate niffler;

use anyhow::{anyhow, Result};
use bio::io::{fasta, fastq};

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
pub fn read_file(path: &Path) -> Result<(Box<dyn io::Read>, niffler::compression::Format)> {
    let raw_in = Box::new(io::BufReader::new(File::open(path)?));

    let (reader, compression) = niffler::get_reader(raw_in).expect("Cannot read input fasta file");

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

/// get_file_type function ---------------------------------------------------

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
pub fn write_to_fa<'a>(
    filename: &'a str,
    compression: niffler::compression::Format,
    out: &'a str,
    record: &'a fasta::Record,
) -> Result<()> {
    let mut file_path = PathBuf::from("");
    file_path.push(out);
    file_path.push(filename);

    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_path)?;
    if compression == niffler::compression::Format::Gzip {
        let handle = niffler::get_writer(
            Box::new(file),
            compression,
            niffler::compression::Level::One,
        )?;
        let mut writer = fasta::Writer::new(handle);
        let _write_res = writer
            .write_record(&record)
            .expect("Cannot write to fasta file");
    } else {
        let handle = io::BufWriter::new(file);
        let mut writer = fasta::Writer::new(handle);
        let _write_res = writer
            .write_record(&record)
            .expect("Cannot write to fasta file");
    }

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
    filename: &'a str,
    compression: niffler::compression::Format,
    out: &'a str,
    record: &'a fastq::Record,
) -> Result<()> {
    let mut file_path = PathBuf::from("");
    file_path.push(out);
    file_path.push(filename);

    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_path)?;

    if compression == niffler::compression::Format::Gzip {
        let handle = niffler::get_writer(
            Box::new(file),
            compression,
            niffler::compression::Level::One,
        )?;
        let mut writer = fastq::Writer::new(handle);
        let _write_res = writer
            .write_record(&record)
            .expect("Cannot write to fasta file");
    } else {
        let handle = io::BufWriter::new(file);
        let mut writer = fastq::Writer::new(handle);
        let _write_res = writer
            .write_record(&record)
            .expect("Cannot write to fasta file");
    }

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
pub fn bc_cmp(bc: &str, seq: &str, mismatch: i32) -> bool {
    let s = &seq[..bc.len()];

    // This wonderful line below compute the number of
    // character mismatch between two strings
    let nb_mismatch: i32 = bc
        .as_bytes()
        .iter()
        .zip(s.as_bytes().iter())
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
pub fn se_fa_demux(
    records: &mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &Barcode,
    mismatch: i32,
    out: &str,
) -> Result<()> {
    let mut ext = "";
    // It is not an obligation to check ext.is_empty in this if-else
    // but it is done here and in subsequent fns to avoid rust warning
    // for used but unread var before assignment
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = "gz";
    } else {
        ext = "";
    }

    let mut nb_records: HashMap<&str, i32> = HashMap::new();
    while let Some(Ok(record)) = records.next() {
        let mut unk = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(record.seq()), mismatch);
            match res {
                true => {
                    if nb_records.get(value[0]) == None {
                        nb_records.insert(value[0], 1);
                    } else {
                        nb_records.insert(value[0], nb_records[value[0]] + 1);
                    }
                    unk = false;
                    write_to_fa(
                        format!("{}.{}", value[0], ext).as_str(),
                        compression,
                        out,
                        &record,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk = true;
                    continue;
                }
            }
        }

        if unk {
            if nb_records.get("unknown.fa") == None {
                nb_records.insert("unknown.fa", 1);
            } else {
                nb_records.insert("unknown.fa", nb_records["unknown.fa"] + 1);
            }
            write_to_fa(
                format!("{}.{}", "unknown.fa", ext).as_str(),
                compression,
                out,
                &record,
            )
            .expect("Cannot write to unknown file");
        }
    }

    let mut sorted: Vec<_> = nb_records.iter().collect();
    sorted.sort_by_key(|a| a.0);

    for (key, value) in sorted.iter() {
        writeln!(io::stdout(), "[INFO] {} contains {} records", key, value)
            .expect("Cannot write to stdout");
    }
    Ok(())
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
pub fn se_fq_demux(
    records: &mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &Barcode,
    mismatch: i32,
    out: &str,
) -> Result<()> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = "gz";
    } else {
        ext = "";
    }

    let mut nb_records: HashMap<&str, i32> = HashMap::new();
    while let Some(Ok(record)) = records.next() {
        let mut unk = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(record.seq()), mismatch);
            match res {
                true => {
                    unk = false;
                    if nb_records.get(value[0]) == None {
                        nb_records.insert(value[0], 1);
                    } else {
                        nb_records.insert(value[0], nb_records[value[0]] + 1);
                    }
                    write_to_fq(
                        format!("{}.{}", value[0], ext).as_str(),
                        compression,
                        out,
                        &record,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk = true;
                    continue;
                }
            }
        }

        if unk {
            if nb_records.get("unknown.fq") == None {
                nb_records.insert("unknown.fq", 1);
            } else {
                nb_records.insert("unknown.fq", nb_records["unknown.fq"] + 1);
            }
            write_to_fq(
                format!("{}.{}", "unknown.fq", ext).as_str(),
                compression,
                out,
                &record,
            )
            .expect("Cannot write to unknown file");
        }
    }
    let mut sorted: Vec<_> = nb_records.iter().collect();
    sorted.sort_by_key(|a| a.0);

    for (key, value) in sorted.iter() {
        writeln!(io::stdout(), "[INFO] {} contains {} records", key, value)
            .expect("Cannot write to stdout");
    }

    Ok(())
}

// pe_fa_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fa_demux(
    forward_records: &mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    reverse_records: &mut fasta::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &Barcode,
    mismatch: i32,
    out: &str,
) -> Result<()> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = "gz";
    } else {
        ext = "";
    }

    let mut nb_records: HashMap<&str, i32> = HashMap::new();
    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut unk = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(f_rec.seq()), mismatch);
            match res {
                true => {
                    unk = false;
                    if nb_records.get(value[0]) == None {
                        nb_records.insert(value[0], 1);
                    } else {
                        nb_records.insert(value[0], nb_records[value[0]] + 1);
                    }
                    write_to_fa(
                        format!("{}.{}", value[0], ext).as_str(),
                        compression,
                        out,
                        &f_rec,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk = true;
                    continue;
                }
            }
        }

        if unk {
            if nb_records.get("unknown_R1.fa") == None {
                nb_records.insert("unknown_R1.fa", 1);
            } else {
                nb_records.insert("unknown_R1.fa", nb_records["unknown_R1.fa"] + 1);
            }
            write_to_fa(
                format!("{}.{}", "unknown_R1.fa", ext).as_str(),
                compression,
                out,
                &f_rec,
            )
            .expect("Cannot write to unknown file");
        }
    }

    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut unk1 = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(r_rec.seq()), mismatch);
            match res {
                true => {
                    if nb_records.get(value[1]) == None {
                        nb_records.insert(value[1], 1);
                    } else {
                        nb_records.insert(value[1], nb_records[value[1]] + 1);
                    }
                    unk1 = false;
                    write_to_fa(
                        format!("{}.{}", value[1], ext).as_str(),
                        compression,
                        out,
                        &r_rec,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk1 = true;
                    continue;
                }
            }
        }

        if unk1 {
            if nb_records.get("unknown_R2.fa") == None {
                nb_records.insert("unknown_R2.fa", 1);
            } else {
                nb_records.insert("unknown_R2.fa", nb_records["unknown_R2.fa"] + 1);
            }
            write_to_fa(
                format!("{}.{}", "unknown_R2.fa", ext).as_str(),
                compression,
                out,
                &r_rec,
            )
            .expect("Cannot write to unknown file");
        }
    }
    let mut sorted: Vec<_> = nb_records.iter().collect();
    sorted.sort_by_key(|a| a.0);

    for (key, value) in sorted.iter() {
        writeln!(io::stdout(), "[INFO] {} contains {} records", key, value)
            .expect("Cannot write to stdout");
    }

    Ok(())
}

// pe_fq_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fq_demux(
    forward_records: &mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    reverse_records: &mut fastq::Records<std::boxed::Box<dyn std::io::Read>>,
    compression: niffler::compression::Format,
    barcode_data: &Barcode,
    mismatch: i32,
    out: &str,
) -> Result<()> {
    let mut ext = "";
    if compression == niffler::compression::Format::Gzip && ext.is_empty() {
        ext = "gz";
    } else {
        ext = "";
    }

    let mut nb_records: HashMap<&str, i32> = HashMap::new();

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut unk = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(f_rec.seq()), mismatch);
            match res {
                true => {
                    unk = false;
                    if nb_records.get(value[0]) == None {
                        nb_records.insert(value[0], 1);
                    } else {
                        nb_records.insert(value[0], nb_records[value[0]] + 1);
                    }
                    write_to_fq(
                        format!("{}.{}", value[0], ext).as_str(),
                        compression,
                        out,
                        &f_rec,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk = true;
                    continue;
                }
            }
        }

        if unk {
            if nb_records.get("unknown_R1.fq") == None {
                nb_records.insert("unknown_R1.fq", 1);
            } else {
                nb_records.insert("unknown_R1.fq", nb_records["unknown_R1.fq"] + 1);
            }
            write_to_fq(
                format!("{}.{}", "unknown_R1.fq", ext).as_str(),
                compression,
                out,
                &f_rec,
            )
            .expect("Cannot write to unknown file");
        }
    }
    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut unk1 = true;
        for (key, value) in barcode_data.iter() {
            let res = bc_cmp(key, &String::from_utf8_lossy(r_rec.seq()), mismatch);
            match res {
                true => {
                    if nb_records.get(value[1]) == None {
                        nb_records.insert(value[1], 1);
                    } else {
                        nb_records.insert(value[1], nb_records[value[1]] + 1);
                    }
                    unk1 = false;
                    write_to_fq(
                        format!("{}.{}", value[1], ext).as_str(),
                        compression,
                        out,
                        &r_rec,
                    )
                    .expect("Cannot write to output file");
                    break;
                }
                false => {
                    unk1 = true;
                    continue;
                }
            }
        }

        if unk1 {
            if nb_records.get("unknown_R2.fq") == None {
                nb_records.insert("unknown_R2.fq", 1);
            } else {
                nb_records.insert("unknown_R2.fq", nb_records["unknown_R2.fq"] + 1);
            }
            write_to_fq(
                format!("{}.{}", "unknown_R2.fq", ext).as_str(),
                compression,
                out,
                &r_rec,
            )
            .expect("Cannot write to unknown file");
        }
    }

    let mut sorted: Vec<_> = nb_records.iter().collect();
    sorted.sort_by_key(|a| a.0);

    for (key, value) in sorted.iter() {
        writeln!(io::stdout(), "[INFO] {} contains {} records", key, value)
            .expect("Cannot write to stdout");
    }

    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::prelude::*;

    // read_file tests ------------------------------------------------------
    #[test]
    fn test_read_file() {
        let (mut reader, compression) = read_file(Path::new("tests/test.fa.gz")).unwrap();
        let mut contents = String::new();
        reader
            .read_to_string(&mut contents)
            .expect("Error during file reading");

        assert_eq!(compression, niffler::compression::Format::Gzip);
        assert_eq!(contents, ">seqID1 desc\nATCGATCGATCGATC\n");
    }

    #[test]
    fn test_read_file_content_is_ok() {
        let (reader, _compression) = read_file(Path::new("tests/test.fa.gz")).unwrap();
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
        let seq = "ATCGATCGATCG";
        let bc = "ATCG";

        assert!(bc_cmp(bc, seq, 0));
    }

    #[test]
    fn test_bc_cmp_not_ok() {
        let bc = "TGCA";
        let seq = "ATCGATCGATCG";

        assert_eq!(bc_cmp(bc, seq, 0), false);
    }

    #[test]
    fn test_bc_cmp_mismatch_ok() {
        let bc = "AACG";
        let seq = "ATCGATCGATCG";

        assert!(bc_cmp(bc, seq, 1));
    }

    #[test]
    fn test_bc_cmp_mismatch_not_ok() {
        let bc = "AACG";
        let seq = "ATCGATCGATCG";

        assert_eq!(bc_cmp(bc, seq, 0), false);
    }

    // se_fa_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fa_demux() {
        let p = Path::new("tests/test2.fa.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("ACCGTA", vec!["id1.fa"]);
        bc_data.insert("ATTGTT", vec!["id2.fa"]);

        assert!(se_fa_demux(&mut records, cmp, &bc_data, 0, out).is_ok());
    }

    #[test]
    fn test_se_fa_demux_m1() {
        let p = Path::new("tests/test2.fa.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("AGCGTA", vec!["id1.fa"]);
        bc_data.insert("ACTGTT", vec!["id2.fa"]);

        assert!(se_fa_demux(&mut records, cmp, &bc_data, 1, out).is_ok());
    }

    #[test]
    fn test_se_fa_demux_m2() {
        let p = Path::new("tests/test2.fa.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fasta::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("GGCGCA", vec!["id1.fa"]);
        bc_data.insert("GCTGCT", vec!["id2.fa"]);

        assert!(se_fa_demux(&mut records, cmp, &bc_data, 2, out).is_ok());
    }

    // se_fq_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fq_demux() {
        let p = Path::new("tests/test2.fq.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("ACCGTA", vec!["id1.fq"]);
        bc_data.insert("ATTGTT", vec!["id2.fq"]);

        assert!(se_fq_demux(&mut records, cmp, &bc_data, 0, out).is_ok());
    }

    #[test]
    fn test_se_fq_demux_m1() {
        let p = Path::new("tests/test2.fq.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("AGCGTA", vec!["id1.fq"]);
        bc_data.insert("AGTGTT", vec!["id2.fq"]);

        assert!(se_fq_demux(&mut records, cmp, &bc_data, 1, out).is_ok());
    }

    #[test]
    fn test_se_fq_demux_m2() {
        let p = Path::new("tests/test2.fq.gz");
        let out = "tests";
        let (fr, cmp) = read_file(&p).expect("Cannot open file");
        let mut records = fastq::Reader::new(fr).records();
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert("AGTGTA", vec!["id1.fq"]);
        bc_data.insert("ACAGTT", vec!["id2.fq"]);

        assert!(se_fq_demux(&mut records, cmp, &bc_data, 2, out).is_ok());
    }

    // write_to_fa tests ----------------------------------------------------
    #[test]
    fn test_write_to_fa_is_ok() {
        let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
        let out = "tests";
        let cmp = niffler::compression::Format::Gzip;
        assert!((write_to_fa("mytmp.fa", cmp, out, &record)).is_ok());

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
        let record = fastq::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG", b"QQQQQQQ");
        let out = "tests";
        let cmp = niffler::compression::Format::Gzip;
        assert!((write_to_fq("mytmp.fq", cmp, out, &record)).is_ok());

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
