// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs::File;
use std::io;
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bio::io::{fasta, fasta::FastaRead, fastq, fastq::FastqRead};

// Get reader and compression format of file
pub fn read_file(
    filename: &str,
) -> Result<(Box<dyn io::Read>, niffler::compression::Format)> {
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).expect("file should be readable"),
    ));

    niffler::get_reader(raw_in).with_context(|| {
        anyhow!("Could not detect compression of file '{}'", filename)
    })
}

pub fn is_fa<P: AsRef<Path> + std::fmt::Debug>(path: &P) -> bool {
    let (stream_in, _) =
        niffler::from_path(path).expect("supplied file should be readable");

    let mut reader = fasta::Reader::new(stream_in);

    let mut record = fasta::Record::new();

    if reader.read(&mut record).is_err() {
        false
    } else {
        record.check().is_ok()
    }
}

pub fn is_fq<P: AsRef<Path> + std::fmt::Debug>(path: &P) -> bool {
    let (stream_in, _) =
        niffler::from_path(path).expect("supplied file should be readable");

    let mut reader = fastq::Reader::new(stream_in);

    let mut record = fastq::Record::new();

    if reader.read(&mut record).is_err() {
        false
    } else {
        record.check().is_ok()
    }
}

pub fn which_format(filename: &str) -> niffler::compression::Format {
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).expect("file should be readable"),
    ));

    let (_, compression) = niffler::sniff(raw_in).expect("cannot");

    compression
}

// FileType structure
#[derive(Debug, PartialEq, Eq)]
pub enum FileType {
    Fasta,
    Fastq,
}

pub fn get_filetype<P: AsRef<Path> + std::fmt::Debug>(path: &P) -> FileType {
    if is_fa(path) {
        FileType::Fasta
    } else {
        FileType::Fastq
    }
}

// Write to provided data to a fasta file in append mode
pub fn write_fa<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fasta::Record,
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)?;

    let mut writer = fasta::Writer::with_capacity(16264, handle);
    writer.write_record(record)?;

    Ok(())
}

// Write to provided data to a fastq file in append mode
pub fn write_fq<'a>(
    file: &'a std::fs::File,
    compression: niffler::compression::Format,
    record: &'a fastq::Record,
    level: niffler::Level,
) -> Result<()> {
    let handle = niffler::get_writer(Box::new(file), compression, level)
        .with_context(|| anyhow!("Could not get file writer"))?;

    let mut writer = fastq::Writer::with_capacity(16264, handle);
    writer
        .write_record(record)
        .expect("file should be wrtiable");

    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::prelude::*;

    #[test]
    fn test_is_fa() {
        let p = Path::new("tests/reads_1.fa.bz2");
        assert!(is_fa(&p));
    }

    #[test]
    fn test_is_fa_error() {
        let p = Path::new("tests/test.fq");
        assert!(!is_fa(&p));
    }

    #[test]
    fn test_is_fq() {
        let p = Path::new("tests/test.fq");
        assert!(is_fq(&p));
    }

    #[test]
    fn test_is_fq_error() {
        let p = Path::new("tests/reads_1.fa");
        assert!(!is_fq(&p));
    }

    #[test]
    fn test_which_format() {
        assert_eq!(
            which_format("tests/test.fa.gz"),
            niffler::compression::Format::Gzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.bz2"),
            niffler::compression::Format::Bzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.xz"),
            niffler::compression::Format::Lzma
        );
        assert_eq!(
            which_format("tests/test.fq.zst"),
            niffler::compression::Format::Zstd
        );
    }

    #[test]
    fn test_get_filetype_fa() {
        let p = Path::new("tests/reads_1.fa");
        assert_eq!(get_filetype(&p), FileType::Fasta);
    }

    #[test]
    fn test_get_filetype_fq() {
        let p = Path::new("tests/test.fq");
        assert_eq!(get_filetype(&p), FileType::Fastq);
    }

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
}
