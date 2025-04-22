// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use anyhow::anyhow;
use fern::colors::ColoredLevelConfig;

use crate::cli;

pub fn setup_logging(quiet: bool) -> anyhow::Result<()> {
    let colors = ColoredLevelConfig::default();
    let level_filter = if quiet {
        log::LevelFilter::Warn
    } else {
        log::LevelFilter::Debug
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

    fern::Dispatch::new()
        .level(level_filter)
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

pub fn create_relpath_from(
    basedir: &Path,
    filename: &str,
    extension: niffler::send::compression::Format,
) -> PathBuf {
    basedir.join(format!("{}{}", filename, to_compression_ext(extension)))
}

// to_niffler_format function
pub fn to_niffler_format(
    format: cli::CompressionFormat,
) -> anyhow::Result<niffler::send::compression::Format> {
    match format {
        cli::CompressionFormat::Gz => Ok(niffler::send::compression::Format::Gzip),
        cli::CompressionFormat::Bz2 => Ok(niffler::send::compression::Format::Bzip),
        cli::CompressionFormat::Xz => Ok(niffler::send::compression::Format::Lzma),
        cli::CompressionFormat::Zst => Ok(niffler::send::compression::Format::Zstd),
    }
}

// Convert niffler compression format to a file extension
pub fn to_compression_ext(compression: niffler::send::compression::Format) -> String {
    match compression {
        niffler::send::compression::Format::Gzip => ".gz".to_string(),
        niffler::send::compression::Format::Bzip => ".bz2".to_string(),
        niffler::send::compression::Format::Lzma => ".xz".to_string(),
        niffler::send::compression::Format::Zstd => ".zst".to_string(),
        niffler::send::compression::Format::No => "".to_string(),
    }
}

// Convert an integer to a niffler::Level
pub fn to_niffler_level(int_level: u8) -> niffler::Level {
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

// Split a &str at each \t
pub fn split_by_tab(string: &str) -> anyhow::Result<Vec<Vec<&str>>> {
    if string.contains('\t') {
        Ok(string
            .lines()
            .map(|line| line.split('\t').collect())
            .collect())
    } else {
        Err(anyhow!("string is not tab-delimited"))
    }
}

// Compare provided barcode with a sequence
pub fn bc_cmp(bc: &[u8], seq: &[u8], mismatch: u8) -> bool {
    // This wonderful line below compute the number of
    // character mismatch between two strings
    bc.iter()
        .zip(seq.iter())
        .map(|(a, b)| (a != b) as u8)
        .sum::<u8>()
        <= mismatch
}

pub fn which_format(filename: &str) -> niffler::send::compression::Format {
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).expect("file should be readable"),
    ));

    let (_, compression) = niffler::send::sniff(raw_in).expect("cannot");

    compression
}

// Write to provided data to a fasta file in append mode
pub fn write_seqs<'a>(
    file: &'a std::fs::File,
    compression: niffler::send::compression::Format,
    record: &'a needletail::parser::SequenceRecord,
    level: niffler::Level,
) -> anyhow::Result<()> {
    let mut handle = niffler::send::get_writer(Box::new(file), compression, level)?;

    match record.format() {
        needletail::parser::Format::Fasta => needletail::parser::write_fasta(
            record.id(),
            &record.seq(),
            &mut handle,
            needletail::parser::LineEnding::Unix,
        )?,
        needletail::parser::Format::Fastq => needletail::parser::write_fastq(
            record.id(),
            &record.seq(),
            record.qual(),
            &mut handle,
            needletail::parser::LineEnding::Unix,
        )?,
    }

    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_relpath_from() {
        assert_eq!(
            create_relpath_from(
                &mut PathBuf::from("path"),
                "file",
                niffler::send::compression::Format::Gzip
            ),
            PathBuf::from("path/file.gz")
        );
    }

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

        assert!(!bc_cmp(bc, seq, 0));
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

        assert!(!bc_cmp(bc, seq, 0));
    }

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
    fn test_split_by_tab_not_ok() {
        let mystring = "HelloWorldEarth\nBrianwasthere";
        assert_eq!(split_by_tab(mystring).is_err(), true);
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
            to_niffler_format(cli::CompressionFormat::Gz).unwrap(),
            niffler::send::compression::Format::Gzip
        );
        assert_eq!(
            to_niffler_format(cli::CompressionFormat::Xz).unwrap(),
            niffler::send::compression::Format::Lzma
        );
        assert_eq!(
            to_niffler_format(cli::CompressionFormat::Bz2).unwrap(),
            niffler::send::compression::Format::Bzip
        );
        assert_eq!(
            to_niffler_format(cli::CompressionFormat::Zst).unwrap(),
            niffler::send::compression::Format::Zstd
        );
    }

    #[test]
    fn test_to_compression_ext() {
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Gzip),
            *".gz"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Lzma),
            *".xz"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Bzip),
            *".bz2"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Zstd),
            *".zst"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::No),
            *""
        );
    }

    #[test]
    fn test_which_format() {
        assert_eq!(
            which_format("tests/test.fa.gz"),
            niffler::send::compression::Format::Gzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.bz2"),
            niffler::send::compression::Format::Bzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.xz"),
            niffler::send::compression::Format::Lzma
        );
        assert_eq!(
            which_format("tests/test.fq.zst"),
            niffler::send::compression::Format::Zstd
        );
    }

    /*
    #[test]
    fn test_write_to_fa_is_ok() {
        let data = b">id_str desc\nATCGCCG";
        let mut reader = noodles::fasta::Reader::new(&data[..]);
        let record = reader.records().next().transpose().unwrap().unwrap();
        let cmp = niffler::compression::Format::Gzip;
        let file = tempfile::tempfile().expect("Cannot create temp file");

        assert!((write_fa(&file, cmp, &record, niffler::Level::One)).is_ok());

        let mut tmpfile =
            tempfile::NamedTempFile::new().expect("Cannot create temp file");
        writeln!(tmpfile, ">id_str desc\nATCGCCG")
            .expect("Cannot write to tmp file");

        let mut fa_records = File::open(tmpfile).map(BufReader::new).map(noodles::fasta::Reader::new).expect("Cannot read file");

        while let Some(Ok(rec)) = fa_records.records().next() {
            assert_eq!(rec.definition().name(), b"id_str");
            assert_eq!(rec.definition().description().unwrap(), b"desc");
            assert_eq!(rec.sequence().as_ref(), b"ATCGCCG");
        }
    }*/
}
