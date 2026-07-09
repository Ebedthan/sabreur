// Copyright 2021-2026 Anicet Ebou.
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
pub fn to_compression_ext(compression: niffler::send::compression::Format) -> &'static str {
    match compression {
        niffler::send::compression::Format::Gzip => ".gz",
        niffler::send::compression::Format::Bzip => ".bz2",
        niffler::send::compression::Format::Lzma => ".xz",
        niffler::send::compression::Format::Zstd => ".zst",
        niffler::send::compression::Format::No => "",
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
        Err(anyhow!("Input string is not tab-delimited"))
    }
}

// Apply a best-match algorithm to find the barcode with the fewest mismatches
pub fn best_barcode_match<'a>(
    barcodes: &[&'a [u8]],
    seq: &[u8],
    max_mismatch: u8,
) -> Option<&'a [u8]> {
    let mut best: Option<(u8, &'a [u8])> = None;
    let mut tied = false;

    for &bc in barcodes {
        let mismatches = bc
            .iter()
            .zip(seq.iter())
            .map(|(a, b)| (a != b) as u8)
            .sum::<u8>();
        if mismatches > max_mismatch {
            continue;
        }

        match best {
            None => best = Some((mismatches, bc)),
            Some((best_count, _)) if mismatches < best_count => {
                best = Some((mismatches, bc));
                tied = false;
            }
            Some((best_count, _)) if mismatches == best_count => tied = true,
            _ => {}
        }
    }

    if tied {
        None
    } else {
        best.map(|(_, bc)| bc)
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

pub fn which_format(filename: &str) -> anyhow::Result<niffler::send::compression::Format> {
    let file = File::open(filename)?;
    let raw_in = Box::new(io::BufReader::new(file));
    let (_, compression) = niffler::send::sniff(raw_in)?;
    Ok(compression)
}

// Write to provided data to a fasta file in append mode
pub fn write_seqs(
    writer: &mut dyn std::io::Write,
    record: &needletail::parser::SequenceRecord,
) -> anyhow::Result<()> {
    match record.format() {
        needletail::parser::Format::Fasta => needletail::parser::write_fasta(
            record.id(),
            record.seq().as_ref(),
            writer,
            needletail::parser::LineEnding::Unix,
        ),
        needletail::parser::Format::Fastq => needletail::parser::write_fastq(
            record.id(),
            record.seq().as_ref(),
            record.qual(),
            writer,
            needletail::parser::LineEnding::Unix,
        ),
    }?;

    Ok(())
}

// Create the persistent, per-output-file writer for a barcode's output
// file: opens the file once and wraps it once with the appropriate
// (de)compressor. Called once per output file at startup, not per record.
//
// Return type is a plain anyhow::Result (error type anyhow::Error), for
// consistency with every other function in this file. niffler::get_writer
// returns Result<_, niffler::Error>; the `?` below converts that into
// anyhow::Error via anyhow's blanket From impl, the same way the `?` on
// `.open(path)` already converts std::io::Error.
pub fn create_output_writer(
    path: &Path,
    compression: niffler::send::compression::Format,
    level: niffler::Level,
) -> anyhow::Result<Box<dyn std::io::Write + Send + 'static>> {
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    Ok(niffler::send::get_writer(
        Box::new(file),
        compression,
        level,
    )?)
}

const IUPAC_CODES: &[u8] = b"ACGTUNRYSWKMBDHV";

// Validates the fields of a barcode file, ensuring it has the correct number of columns,
// contains only valid IUPAC nucleotide codes, and has unique barcodes.
pub fn validate_barcode_fields(fields: &[Vec<&str>], paired: bool) -> anyhow::Result<()> {
    if fields.is_empty() {
        return Err(anyhow!("Barcode file is empty"));
    }

    let expected_cols = if paired { 3 } else { 2 };
    let mut seen_barcodes: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut bc_len: Option<usize> = None;

    for (idx, row) in fields.iter().enumerate() {
        let line_no = idx + 1;

        if row.len() != expected_cols {
            return Err(anyhow!(
                "Barcode file line {}: expected {} tab-separated column(s) ({}), found {}",
                line_no,
                expected_cols,
                if paired {
                    "barcode, forward sample name, reverse sample name"
                } else {
                    "barcode, sample name"
                },
                row.len()
            ));
        }

        let barcode = row[0];

        if row.iter().any(|field| field.is_empty()) {
            return Err(anyhow!(
                "Barcode file line {}: contains an empty field",
                line_no
            ));
        }

        if barcode.eq_ignore_ascii_case("XXX") {
            return Err(anyhow!("Barcode file line {}: 'XXX' is reserved internally for unknown/unmatched reads and cannot be used as a barcode", line_no));
        }

        if !barcode
            .bytes()
            .all(|b| IUPAC_CODES.contains(&b.to_ascii_uppercase()))
        {
            return Err(anyhow!("Barcode file line {}: barcode '{}' contains characters that are not valid IUPAC nucleotide codes", line_no, barcode));
        }

        match bc_len {
            None => bc_len = Some(barcode.len()),
            Some(len) if len != barcode.len() => {
                return Err(anyhow!("Barcode file line {}: barcode '{}' has length {}, but earlier barcodes in this file have length {}. All barcodes must be the same length.",
                    line_no, barcode, barcode.len(), len));
            }
            _ => {}
        }

        if !seen_barcodes.insert(barcode) {
            return Err(anyhow!(
                "Barcode file line {}: barcode '{}' appears more than once",
                line_no,
                barcode
            ));
        }
    }

    Ok(())
}

// Tests -----------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_relpath_from() {
        assert_eq!(
            create_relpath_from(
                &PathBuf::from("path"),
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
        assert!(split_by_tab(mystring).is_err());
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
            ".gz"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Lzma),
            ".xz"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Bzip),
            ".bz2"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::Zstd),
            ".zst"
        );
        assert_eq!(
            to_compression_ext(niffler::send::compression::Format::No),
            ""
        );
    }

    #[test]
    fn test_which_format() {
        assert_eq!(
            which_format("tests/test.fa.gz").unwrap(),
            niffler::send::compression::Format::Gzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.bz2").unwrap(),
            niffler::send::compression::Format::Bzip
        );
        assert_eq!(
            which_format("tests/reads_1.fa.xz").unwrap(),
            niffler::send::compression::Format::Lzma
        );
        assert_eq!(
            which_format("tests/test.fq.zst").unwrap(),
            niffler::send::compression::Format::Zstd
        );
    }

    // -- best_barcode_match ---------------------------------------------

    #[test]
    fn test_best_barcode_match_exact_wins_over_mismatch() {
        // "AACGTA" is 1 mismatch from ACCGTA and 2+ from ATTGTT; ACCGTA
        // is the correct pick even though a naive first-found search
        // (depending on hashmap order) could have picked either.
        let barcodes: Vec<&[u8]> = vec![b"ACCGTA", b"ATTGTT"];
        let seq = b"ACCGTA";
        assert_eq!(
            best_barcode_match(&barcodes, seq, 1),
            Some(b"ACCGTA".as_ref())
        );
    }

    #[test]
    fn test_best_barcode_match_prefers_lower_mismatch_count() {
        // seq is 1 mismatch from "AACGTA" and 3 mismatches from "TTTGTA";
        // the lower mismatch count must win even though both are within
        // threshold.
        let barcodes: Vec<&[u8]> = vec![b"AACGTA", b"TTTGTA"];
        let seq = b"ACCGTA";
        assert_eq!(
            best_barcode_match(&barcodes, seq, 3),
            Some(b"AACGTA".as_ref())
        );
    }

    #[test]
    fn test_best_barcode_match_ambiguous_tie_returns_none() {
        // seq is exactly 1 mismatch from both barcodes, a genuine tie.
        // The correct behavior is None (caller routes to unknown), not an
        // arbitrary pick of whichever came first.
        let barcodes: Vec<&[u8]> = vec![b"AACGTA", b"ACCGTT"];
        let seq = b"ACCGTA";
        assert_eq!(best_barcode_match(&barcodes, seq, 1), None);
    }

    #[test]
    fn test_best_barcode_match_no_candidate_within_threshold() {
        let barcodes: Vec<&[u8]> = vec![b"TTTTTT"];
        let seq = b"ACCGTA";
        assert_eq!(best_barcode_match(&barcodes, seq, 1), None);
    }

    #[test]
    fn test_best_barcode_match_single_exact_candidate() {
        let barcodes: Vec<&[u8]> = vec![b"ACCGTA"];
        let seq = b"ACCGTA";
        assert_eq!(
            best_barcode_match(&barcodes, seq, 0),
            Some(b"ACCGTA".as_ref())
        );
    }

    // -- validate_barcode_fields ------------------------------------------

    #[test]
    fn test_validate_barcode_fields_valid_se() {
        let fields = vec![vec!["ACCGTA", "sample1"], vec!["ATTGTT", "sample2"]];
        assert!(validate_barcode_fields(&fields, false).is_ok());
    }

    #[test]
    fn test_validate_barcode_fields_valid_pe() {
        let fields = vec![
            vec!["ACCGTA", "sample1_R1", "sample1_R2"],
            vec!["ATTGTT", "sample2_R1", "sample2_R2"],
        ];
        assert!(validate_barcode_fields(&fields, true).is_ok());
    }

    #[test]
    fn test_validate_barcode_fields_empty_file() {
        let fields: Vec<Vec<&str>> = vec![];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_wrong_column_count() {
        // PE mode expects 3 columns; this row only has 2.
        let fields = vec![vec!["ACCGTA", "sample1"]];
        assert!(validate_barcode_fields(&fields, true).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_empty_field_rejected() {
        let fields = vec![vec!["ACCGTA", ""]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_xxx_barcode_rejected() {
        let fields = vec![vec!["XXX", "sample1"]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_xxx_rejected_case_insensitive() {
        let fields = vec![vec!["xxx", "sample1"]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_non_iupac_characters_rejected() {
        let fields = vec![vec!["ACCZTA", "sample1"]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_mixed_length_rejected() {
        let fields = vec![vec!["ACCGTA", "sample1"], vec!["ACC", "sample2"]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    #[test]
    fn test_validate_barcode_fields_duplicate_rejected() {
        let fields = vec![vec!["ACCGTA", "sample1"], vec!["ACCGTA", "sample2"]];
        assert!(validate_barcode_fields(&fields, false).is_err());
    }

    // -- create_output_writer ---------------------------------------------

    #[test]
    fn test_create_output_writer_roundtrip_uncompressed() {
        use std::io::Read;

        let tmp = tempfile::NamedTempFile::new().expect("Cannot create temp file");
        let path = tmp.path().to_path_buf();

        {
            let mut writer = create_output_writer(
                &path,
                niffler::send::compression::Format::No,
                niffler::Level::One,
            )
            .expect("create_output_writer failed");
            writer.write_all(b">read1\nACGT\n").expect("write failed");
            writer.flush().expect("flush failed");
        }

        let mut content = String::new();
        std::fs::File::open(&path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();
        assert_eq!(content, ">read1\nACGT\n");
    }

    #[test]
    fn test_create_output_writer_roundtrip_gzip() {
        let tmp = tempfile::NamedTempFile::new().expect("Cannot create temp file");
        let path = tmp.path().to_path_buf();

        {
            let mut writer = create_output_writer(
                &path,
                niffler::send::compression::Format::Gzip,
                niffler::Level::One,
            )
            .expect("create_output_writer failed");
            writer.write_all(b">read1\nACGT\n").expect("write failed");
            // Dropping here finalizes the gzip stream (see the Drop-based
            // finalization caveat discussed for niffler's boxed writers).
        }

        let (_, format) =
            niffler::send::sniff(Box::new(io::BufReader::new(File::open(&path).unwrap())))
                .expect("sniff failed");
        assert_eq!(format, niffler::send::compression::Format::Gzip);
    }
}
