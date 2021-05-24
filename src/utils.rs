// Copyright 2021 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self};

extern crate niffler;

use anyhow::{anyhow, Context, Result};
use bio::io::{fasta, fastq};

use crate::error;

// to_niffler_format function
pub fn to_niffler_format(format: &str) -> Result<niffler::compression::Format> {
    if format == "gz" {
        Ok(niffler::compression::Format::Gzip)
    } else if format == "bz2" {
        Ok(niffler::compression::Format::Bzip)
    } else if format == "xz" {
        Ok(niffler::compression::Format::Lzma)
    } else {
        Ok(niffler::compression::Format::No)
    }
}

// to_compression_ext function

/// Convert niffler compression format to a file extension
///
/// # Example
/// ```rust
/// let compression = niffler::compression::Format::Gzip;
/// let ext = to_compression_ext(compression);
/// assert_eq!(ext, ".gz");
/// ```
///
pub fn to_compression_ext(compression: niffler::compression::Format) -> String {
    match compression {
        niffler::compression::Format::Gzip => ".gz".to_string(),
        niffler::compression::Format::Bzip => ".bz2".to_string(),
        niffler::compression::Format::Lzma => ".xz".to_string(),
        niffler::compression::Format::No => "".to_string(),
    }
}

// to_niffler_level function

/// Convert an integer to a niffler::Level
///
/// # Example
/// ```rust
/// assert_eq!(to_niffler_level(3), niffler::Level::Three);
/// ```
///
pub fn to_niffler_level(int_level: i32) -> niffler::Level {
    match int_level {
        1 => niffler::Level::One,
        2 => niffler::Level::Two,
        3 => niffler::Level::Three,
        4 => niffler::Level::Four,
        5 => niffler::Level::Five,
        6 => niffler::Level::Six,
        7 => niffler::Level::Seven,
        8 => niffler::Level::Eight,
        9 => niffler::Level::Nine,
        _ => niffler::Level::One,
    }
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
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).with_context(|| error::Error::CantReadFile {
            filename: filename.to_string(),
        })?,
    ));

    niffler::get_reader(raw_in).with_context(|| {
        anyhow!("Could not detect compression of file '{}'", filename)
    })
}

/// get_file_type function --------------------------------------------------

// FileType structure
#[derive(Debug, PartialEq)]
pub enum FileType {
    Fasta,
    Fastq,
    None,
}

/// Get file type from filename
///
/// # Example
/// ```rust
/// let filename = "myfile.fq";
/// let file_type = get_file_type(filename);
/// ```
///
pub fn get_file_type(
    filename: &str,
) -> Result<(FileType, niffler::compression::Format)> {
    let mut filext = FileType::None;

    if filext == FileType::None && filename.contains(".fastq")
        || filename.contains(".fq")
    {
        filext = FileType::Fastq;
    } else if filename.contains(".fasta")
        || filename.contains(".fa")
        || filename.contains(".fas")
    {
        filext = FileType::Fasta;
    } else {
        filext = FileType::None;
    }

    let (_reader, comp) = read_file(filename)?;

    Ok((filext, comp))
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
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)?;

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
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)
        .with_context(|| anyhow!("Could not get file writer"))?;

    let mut writer = fastq::Writer::new(handle);
    let _write_res = writer.write_record(&record).with_context(|| {
        error::Error::CantWriteFile {
            filename: "output file".to_string(),
        }
    })?;

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
    forward: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool)> {
    let (forward_reader, mut compression) =
        read_file(forward).with_context(|| error::Error::ReadingError {
            filename: forward.to_string(),
        })?;
    let mut forward_records = fasta::Reader::new(forward_reader).records();

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    let mut is_unk_empty = true;

    if format != niffler::compression::Format::No {
        compression = format;
    }

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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
            }
            None => {
                is_unk_empty = false;
                write_to_fa(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                    compression,
                    &f_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
            }
        }
    }

    Ok((nb_records, is_unk_empty))
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
    forward: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool)> {
    let (forward_reader, mut compression) =
        read_file(forward).with_context(|| error::Error::ReadingError {
            filename: forward.to_string(),
        })?;
    let mut forward_records = fastq::Reader::new(forward_reader).records();

    if format != niffler::compression::Format::No {
        compression = format;
    }

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    let mut is_unk_empty = true;

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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
            }
            None => {
                is_unk_empty = false;
                write_to_fq(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                    compression,
                    &f_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
            }
        }
    }

    Ok((nb_records, is_unk_empty))
}

// pe_fa_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fa_demux<'a>(
    forward: &'a str,
    reverse: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool, bool)> {
    let (forward_reader, mut compression) =
        read_file(forward).with_context(|| error::Error::ReadingError {
            filename: forward.to_string(),
        })?;
    let mut forward_records = fasta::Reader::new(forward_reader).records();

    let (reverse_reader, _compression) =
        read_file(reverse).with_context(|| error::Error::ReadingError {
            filename: reverse.to_string(),
        })?;
    let mut reverse_records = fasta::Reader::new(reverse_reader).records();

    if format != niffler::compression::Format::No {
        compression = format;
    }

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    let mut is_unk_r1_empty = true;
    let mut is_unk_r2_empty = true;

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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
            }
            None => {
                is_unk_r1_empty = false;
                write_to_fa(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                    compression,
                    &f_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
            }
            None => {
                is_unk_r2_empty = false;
                write_to_fa(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[1],
                    compression,
                    &r_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fasta,
                    }
                })?;
            }
        }
    }

    Ok((nb_records, is_unk_r1_empty, is_unk_r2_empty))
}

// pe_fq_demux function -----------------------------------------------------

/// Demultiplex a fasta::Record of paired-end file
///
///
///
pub fn pe_fq_demux<'a>(
    forward: &'a str,
    reverse: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool, bool)> {
    let (forward_reader, mut compression) =
        read_file(forward).with_context(|| error::Error::ReadingError {
            filename: forward.to_string(),
        })?;
    let mut forward_records = fastq::Reader::new(forward_reader).records();

    let (reverse_reader, _compression) =
        read_file(reverse).with_context(|| error::Error::ReadingError {
            filename: reverse.to_string(),
        })?;
    let mut reverse_records = fastq::Reader::new(reverse_reader).records();

    if format != niffler::compression::Format::No {
        compression = format;
    }

    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();
    let bc_len = my_vec[0].len();

    let mut is_unk_r1_empty = true;
    let mut is_unk_r2_empty = true;

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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
            }
            None => {
                is_unk_r1_empty = false;
                write_to_fq(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                    compression,
                    &f_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
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
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
            }
            None => {
                is_unk_r2_empty = false;
                write_to_fq(
                    &barcode_data.get(&"XXX".as_bytes()).unwrap()[1],
                    compression,
                    &r_rec,
                    level,
                )
                .with_context(|| {
                    error::Error::WritingErrorNoFilename {
                        format: FileType::Fastq,
                    }
                })?;
            }
        }
    }

    Ok((nb_records, is_unk_r1_empty, is_unk_r2_empty))
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::OpenOptions;
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fa_demux(
            "tests/test2.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m1() {
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fa_demux(
            "tests/test2.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m2() {
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fa_demux(
            "tests/test2.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
            &mut nb_records
        )
        .is_ok());
    }

    // se_fq_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fq_demux() {
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fq_demux(
            "tests/test2.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m1() {
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fq_demux(
            "tests/test2.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m2() {
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
        let unknown_file1 = OpenOptions::new()
            .create(true)
            .append(true)
            .open("tests/unk1.fa")
            .expect("cannot open file");
        bc_data.insert(b"ACCGTA", vec![file1]);
        bc_data.insert(b"ATTGTT", vec![file2]);
        bc_data.insert(b"XXX", vec![unknown_file1]);

        assert!(se_fq_demux(
            "tests/test2.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
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
        assert!((write_to_fa(&file, cmp, &record, niffler::Level::One)).is_ok());

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
            .open("tests/mytmp.fq")
            .expect("cannot open file");
        assert!((write_to_fq(&file, cmp, &record, niffler::Level::One)).is_ok());

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
        let f1 = "tests/reads_1.fa";
        let f3 = "tests/reads_1.fa.gz";
        let f5 = "tests/reads_1.fa.xz";

        assert_eq!(
            get_file_type(f1).unwrap(),
            (FileType::Fasta, niffler::compression::Format::No)
        );
        assert_eq!(
            get_file_type(f3).unwrap(),
            (FileType::Fasta, niffler::compression::Format::Gzip)
        );
        assert_eq!(
            get_file_type(f5).unwrap(),
            (FileType::Fasta, niffler::compression::Format::Lzma)
        );
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
