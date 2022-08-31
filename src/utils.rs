// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::PathBuf;

extern crate chrono;
extern crate fern;
extern crate niffler;

use crate::error;
use anyhow::{anyhow, Context, Result};
use bio::io::{fasta, fastq};
use fern::colors::ColoredLevelConfig;

pub fn setup_logging(quiet: bool) -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::default();
    let mut base_config = fern::Dispatch::new();

    base_config = match quiet {
        // if user required quietness let only output warning messages
        // or messages more severe than warnings
        true => base_config.level(log::LevelFilter::Warn),
        // if quietness is not specified which implies verbosity is allowed
        // output
        false => base_config.level(log::LevelFilter::Debug),
    };

    // Separate file config so we can include year, month and day in file logs
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
        .chain(fern::log_file("sabreur.log")?);

    let stdout_config = fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{}][{}] {}",
                chrono::Local::now().format("%H:%M:%S"),
                colors.color(record.level()),
                message
            ))
        })
        .chain(io::stdout());

    base_config
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

pub fn create_relpath_from(input: Vec<&str>) -> Result<PathBuf> {
    let mut path = PathBuf::from("");

    for element in input {
        path.push(element);
    }

    Ok(path)
}

// to_niffler_format function
pub fn to_niffler_format(format: &str) -> Result<niffler::compression::Format> {
    match format {
        "gz" => Ok(niffler::compression::Format::Gzip),
        "bz2" => Ok(niffler::compression::Format::Bzip),
        "xz" => Ok(niffler::compression::Format::Lzma),
        "zst" => Ok(niffler::compression::Format::Zstd),
        _ => Ok(niffler::compression::Format::No),
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
        niffler::compression::Format::Zstd => ".zst".to_string(),
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
pub fn write_fa<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fasta::Record,
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)?;

    let mut writer = fasta::Writer::with_capacity(16264, handle);
    let _write_res = writer.write_record(record)?;

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
pub fn write_fq<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fastq::Record,
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)
        .with_context(|| anyhow!("Could not get file writer"))?;

    let mut writer = fastq::Writer::with_capacity(16264, handle);
    let _write_res = writer.write_record(record).with_context(|| {
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
    file: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool)> {
    // Get fasta file reader and compression mode
    let (reader, mut compression) = read_file(file)?;

    // Get records
    let mut records = fasta::Reader::new(reader).records();

    // Clone barcode values in barcode_data structure for future iteration
    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();

    // Get barcode length
    let bc_len = my_vec[0].len();

    // Initialize unknown file as empty
    let mut is_unk_empty = true;

    // Change output compression format to user wanted compression
    // format if specified by --format option
    if format != niffler::compression::Format::No {
        compression = format;
    }

    while let Some(Ok(record)) = records.next() {
        // Match sequence and barcode with mismatch
        // and return matched barcode. We first use
        // let iter = my_vec.iter() to further stop
        // the find at first match.
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &record.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fa(
                &barcode_data.get(i).unwrap()[0],
                compression,
                &record,
                level,
            )
            .with_context(|| {
                error::Error::WritingErrorNoFilename {
                    format: FileType::Fasta,
                }
            })?;
        } else {
            is_unk_empty = false;
            write_fa(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &record,
                level,
            )
            .with_context(|| {
                error::Error::WritingErrorNoFilename {
                    format: FileType::Fasta,
                }
            })?;
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
    file: &'a str,
    format: niffler::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: i32,
    nb_records: &'a mut HashMap<&'a [u8], i32>,
) -> Result<(&'a mut HashMap<&'a [u8], i32>, bool)> {
    // Get fastq file reader and compression mode
    let (reader, mut compression) = read_file(file)?;

    // Get records
    let mut records = fastq::Reader::new(reader).records();

    // Clone barcode values in barcode_data structure for future iteration
    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();

    // Get barcode length
    let bc_len = my_vec[0].len();

    // Initialize unknown file as empty
    let mut is_unk_empty = true;

    // Change output compression format to user wanted compression
    // format if specified by --format option
    if format != niffler::compression::Format::No {
        compression = format;
    }

    while let Some(Ok(record)) = records.next() {
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &record.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fq(
                &barcode_data.get(i).unwrap()[0],
                compression,
                &record,
                level,
            )
            .with_context(|| {
                error::Error::WritingErrorNoFilename {
                    format: FileType::Fastq,
                }
            })?;
        } else {
            is_unk_empty = false;
            write_fq(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &record,
                level,
            )
            .with_context(|| {
                error::Error::WritingErrorNoFilename {
                    format: FileType::Fastq,
                }
            })?;
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
) -> Result<(&'a mut HashMap<&'a [u8], i32>, String)> {
    // Get fasta files reader and compression modes
    let (forward_reader, mut compression) = read_file(forward)?;

    let (reverse_reader, _compression) = read_file(reverse)?;

    // Get records
    let mut forward_records = fasta::Reader::new(forward_reader).records();
    let mut reverse_records = fasta::Reader::new(reverse_reader).records();

    // Clone barcode values in barcode_data structure for future iteration
    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();

    // Get barcode length
    let bc_len = my_vec[0].len();

    // Initialize unknown files as empty
    let mut unk1_empty = "true";
    let mut unk2_empty = "true";

    // Change output compression format to user wanted compression
    // format if specified by --format option
    if format != niffler::compression::Format::No {
        compression = format;
    }

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fa(
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
        } else {
            unk1_empty = "false";
            write_fa(
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

    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &r_rec.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fa(
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
        } else {
            unk2_empty = "false";
            write_fa(
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
    let final_str = format!("{}{}", unk1_empty, unk2_empty);
    Ok((nb_records, final_str))
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
) -> Result<(&'a mut HashMap<&'a [u8], i32>, String)> {
    // Get fasta files reader and compression modes
    let (forward_reader, mut compression) = read_file(forward)?;

    let (reverse_reader, _compression) = read_file(reverse)?;

    // Get records
    let mut forward_records = fastq::Reader::new(forward_reader).records();
    let mut reverse_records = fastq::Reader::new(reverse_reader).records();

    // Clone barcode values in barcode_data structure for future iteration
    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();

    // Get barcode length
    let bc_len = my_vec[0].len();

    // Initialize unknown files as empty
    let mut unk1_empty = "true";
    let mut unk2_empty = "true";

    // Change output compression format to user wanted compression
    // format if specified by --format option
    if format != niffler::compression::Format::No {
        compression = format;
    }

    while let Some(Ok(f_rec)) = forward_records.next() {
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &f_rec.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fq(
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
        } else {
            unk1_empty = "false";
            write_fq(
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

    while let Some(Ok(r_rec)) = reverse_records.next() {
        let mut iter = my_vec.iter();
        let matched_barcode =
            iter.find(|&&x| bc_cmp(x, &r_rec.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_fq(
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
        } else {
            unk2_empty = "false";
            write_fq(
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
    let final_str = format!("{}{}", unk1_empty, unk2_empty);
    Ok((nb_records, final_str))
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::prelude::*;

    // create_relpath_from --------------------------------------------------
    #[test]
    fn test_create_relpath_from() {
        assert_eq!(
            create_relpath_from(["path", "to", "file"].to_vec()).unwrap(),
            PathBuf::from("path/to/file")
        );
    }

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

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fa_demux(
            "tests/test.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_trim() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fa_demux(
            "tests/test.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fa_demux(
            "tests/test.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_fa_demux_m2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fa_demux(
            "tests/test.fa.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
            &mut nb_records,
        )
        .is_ok());
    }

    // se_fq_demux tests ----------------------------------------------------
    #[test]
    fn test_se_fq_demux() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fq_demux(
            "tests/test.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fq_demux(
            "tests/test.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_fq_demux_m2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], i32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_fq_demux(
            "tests/test.fq.gz",
            niffler::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
            &mut nb_records,
        )
        .is_ok());
    }

    // write_to_fa tests ----------------------------------------------------
    #[test]
    fn test_write_to_fa_is_ok() {
        let record =
            fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
        let cmp = niffler::compression::Format::Gzip;
        let file = tempfile::tempfile().expect("Cannot create temp file");

        assert!((write_fa(&file, cmp, &record, niffler::Level::One)).is_ok());

        let mut tmpfile =
            tempfile::NamedTempFile::new().expect("Cannot create temp file");
        writeln!(tmpfile, ">id_str desc\nATCGCCG")
            .expect("Cannot write to tmp file");

        let mut fa_records = fasta::Reader::from_file(tmpfile)
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
        let file = tempfile::tempfile().expect("Cannot create temp file");

        assert!((write_fq(&file, cmp, &record, niffler::Level::One)).is_ok());

        let mut tmpfile =
            tempfile::NamedTempFile::new().expect("Cannot create temp file");
        writeln!(tmpfile, ">id_str desc\nATCGCCG\n+\nQQQQQQQ")
            .expect("Cannot write to tmp file");

        let mut fa_records = fastq::Reader::from_file(tmpfile)
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

    #[test]
    fn test_to_niffler_level() {
        assert_eq!(to_niffler_level(1), niffler::Level::One);
        assert_eq!(to_niffler_level(2), niffler::Level::Two);
        assert_eq!(to_niffler_level(3), niffler::Level::Three);
        assert_eq!(to_niffler_level(4), niffler::Level::Four);
        assert_eq!(to_niffler_level(5), niffler::Level::Five);
        assert_eq!(to_niffler_level(6), niffler::Level::Six);
        assert_eq!(to_niffler_level(7), niffler::Level::Seven);
        assert_eq!(to_niffler_level(8), niffler::Level::Eight);
        assert_eq!(to_niffler_level(9), niffler::Level::Nine);
    }

    #[test]
    fn test_to_niffler_format() {
        assert_eq!(
            to_niffler_format("gz").unwrap(),
            niffler::compression::Format::Gzip
        );
        assert_eq!(
            to_niffler_format("xz").unwrap(),
            niffler::compression::Format::Lzma
        );
        assert_eq!(
            to_niffler_format("bz2").unwrap(),
            niffler::compression::Format::Bzip
        );
        assert_eq!(
            to_niffler_format("txt").unwrap(),
            niffler::compression::Format::No
        );
    }

    #[test]
    fn test_to_compression_ext() {
        assert_eq!(
            to_compression_ext(niffler::compression::Format::Gzip),
            *".gz"
        );
        assert_eq!(
            to_compression_ext(niffler::compression::Format::Lzma),
            *".xz"
        );
        assert_eq!(
            to_compression_ext(niffler::compression::Format::Bzip),
            *".bz2"
        );
        assert_eq!(to_compression_ext(niffler::compression::Format::No), *"");
    }
}
