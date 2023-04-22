// Copyright 2021-2023 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::io;
use std::path::PathBuf;

extern crate chrono;
extern crate fern;
extern crate niffler;

use anyhow::{anyhow, Result};
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

pub fn create_relpath_from(input: Vec<&str>) -> PathBuf {
    let mut path = PathBuf::from("");

    input.iter().for_each(|x| path.push(x));

    path
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

// Convert niffler compression format to a file extension
pub fn to_compression_ext(compression: niffler::compression::Format) -> String {
    match compression {
        niffler::compression::Format::Gzip => ".gz".to_string(),
        niffler::compression::Format::Bzip => ".bz2".to_string(),
        niffler::compression::Format::Lzma => ".xz".to_string(),
        niffler::compression::Format::Zstd => ".zst".to_string(),
        niffler::compression::Format::No => "".to_string(),
    }
}

// Convert an integer to a niffler::Level
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

// Split a &str at each \t
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

// Compare provided barcode with a sequence
pub fn bc_cmp(bc: &[u8], seq: &[u8], mismatch: i32) -> bool {
    // This wonderful line below compute the number of
    // character mismatch between two strings
    bc.iter()
        .zip(seq.iter())
        .map(|(a, b)| (a != b) as i32)
        .sum::<i32>()
        <= mismatch
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_relpath_from() {
        assert_eq!(
            create_relpath_from(["path", "to", "file"].to_vec()),
            PathBuf::from("path/to/file")
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
            to_niffler_format("zst").unwrap(),
            niffler::compression::Format::Zstd
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
        assert_eq!(
            to_compression_ext(niffler::compression::Format::Zstd),
            *".zst"
        );
        assert_eq!(to_compression_ext(niffler::compression::Format::No), *"");
    }
}
