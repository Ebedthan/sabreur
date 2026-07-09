// Copyright 2021-2026 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;

use crate::utils::{best_barcode_match, write_seqs};

pub type Barcode<'a> = HashMap<&'a [u8], Vec<std::fs::File>>;

/// A function to demultiplex a FASTA/FASTQ file
pub fn se_demux<'a>(
    file: &'a str,
    mut format: niffler::send::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode<'a>,
    mismatch: u8,
    nb_records: &'a mut HashMap<&'a [u8], u32>,
) -> anyhow::Result<(&'a mut HashMap<&'a [u8], u32>, bool)> {
    // Prepare decompression stream
    let (reader, original_format) = niffler::send::from_path(file)?;
    let mut reader = needletail::parse_fastx_reader(reader)?;

    // Use user-specified compression format if set
    if format == niffler::send::compression::Format::No {
        format = original_format;
    }

    // Cache barcode keys (avoid repeated hashmap lookups)
    let mut barcodes: Vec<&[u8]> = barcode_data
        .keys()
        .copied()
        .filter(|&k| k != b"XXX".as_ref())
        .collect();
    barcodes.sort_unstable();

    let Some(&first_key) = barcodes.first() else {
        return Err(anyhow::anyhow!("Barcode data is empty"));
    };

    let bc_len = first_key.len();

    // Track whether unknown file has data
    let mut is_unk_empty = true;

    // Get handle for unkwnon barcode file
    let unknown_writer = barcode_data
        .get(b"XXX".as_ref())
        .ok_or_else(|| anyhow::anyhow!("Missing 'XXX' fallback barcode entry in barcode_data"))?[0]
        .try_clone()?;

    // Process each record
    while let Some(record) = reader.next() {
        let record = record?;

        let seq = record.seq();
        let matched = (seq.len() >= bc_len)
            .then(|| best_barcode_match(&barcodes, &seq[..bc_len], mismatch))
            .flatten();

        match matched {
            Some(bc) => {
                *nb_records.entry(bc).or_insert(0) += 1;
                write_seqs(&barcode_data[bc][0], format, &record, level)?;
            }
            None => {
                is_unk_empty = false;
                write_seqs(&unknown_writer, format, &record, level)?;
            }
        }
    }
    Ok((nb_records, is_unk_empty))
}

/// A function to demultiplex a pair of FASTA/FASTQ files
pub fn pe_demux<'a>(
    forward: &'a str,
    reverse: &'a str,
    format: niffler::send::compression::Format,
    level: niffler::Level,
    barcode_data: &'a Barcode,
    mismatch: u8,
    nb_records: &'a mut HashMap<&'a [u8], u32>,
) -> anyhow::Result<(&'a mut HashMap<&'a [u8], u32>, String)> {
    // Get fasta files reader and compression modes
    let (forward_reader, mut compression) = niffler::send::from_path(forward)?;
    let (reverse_reader, _) = niffler::send::from_path(reverse)?;

    // Get records
    let mut forward_fastx_reader = needletail::parse_fastx_reader(forward_reader)?;
    let mut reverse_fastx_reader = needletail::parse_fastx_reader(reverse_reader)?;

    // Get barcode information once
    // pe_demux
    let mut barcodes: Vec<&[u8]> = barcode_data
        .keys()
        .copied()
        .filter(|&k| k != b"XXX".as_ref())
        .collect();
    barcodes.sort_unstable();

    let Some(&first_barcode) = barcodes.first() else {
        return Err(anyhow::anyhow!("Barcode data is empty"));
    };

    let bc_len = first_barcode.len();
    let unknown_key = "XXX".as_bytes();
    let unknown_files = barcode_data.get(unknown_key).unwrap();

    // Change output compression format if specified
    if format != niffler::send::compression::Format::No {
        compression = format;
    }

    // Process forward and reverse reads together. Using a lockstep so that mates
    // always stay paired. The barcode is looked up from the forward read only (R1),
    // and the resulting alignment (a barcode, or "unknown") is applied to BOTH mates.
    // This keeps output files synchronized.
    let mut unk1_empty = true;
    let mut unk2_empty = true;

    loop {
        let fwd_item = forward_fastx_reader.next();
        let rev_item = reverse_fastx_reader.next();

        let (fwd_record, rev_record) = match (fwd_item, rev_item) {
            (Some(f), Some(r)) => (f?, r?),
            (None, None) => break,
            _ => {
                return Err(anyhow::anyhow!(
                    "Forward and reverse files have a different number of \
                    records. They are not properly paired."
                ));
            }
        };

        let fwd_seq = fwd_record.seq();
        let matched = (fwd_seq.len() >= bc_len)
            .then(|| best_barcode_match(&barcodes, &fwd_seq[..bc_len], mismatch))
            .flatten();

        match matched {
            Some(bc) => {
                *nb_records.entry(bc).or_insert(0) += 1;
                write_seqs(&barcode_data[bc][0], compression, &fwd_record, level)?;
                write_seqs(&barcode_data[bc][1], compression, &rev_record, level)?;
            }
            None => {
                unk1_empty = false;
                unk2_empty = false;
                write_seqs(&unknown_files[0], compression, &fwd_record, level)?;
                write_seqs(&unknown_files[1], compression, &rev_record, level)?;
            }
        }
    }

    // Create result string more efficiently
    let final_str = format!("{unk1_empty}{unk2_empty}");
    Ok((nb_records, final_str))
}

// Tests ----------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_se_demux_1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fa.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_trim() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fa.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_m1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fa.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_m2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fa.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fq.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            0,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_m3() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fq.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            1,
            &mut nb_records,
        )
        .is_ok());
    }

    #[test]
    fn test_se_demux_m4() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        let forward = tempfile::tempfile().expect("Cannot create temp file");
        let reverse = tempfile::tempfile().expect("Cannot create temp file");
        let unknown = tempfile::tempfile().expect("Cannot create temp file");

        bc_data.insert(b"ACCGTA", vec![forward]);
        bc_data.insert(b"ATTGTT", vec![reverse]);
        bc_data.insert(b"XXX", vec![unknown]);

        assert!(se_demux(
            "tests/test.fq.gz",
            niffler::send::compression::Format::Gzip,
            niffler::Level::One,
            &bc_data,
            2,
            &mut nb_records,
        )
        .is_ok());
    }
}
