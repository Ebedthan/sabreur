// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;

use crate::utils::{bc_cmp, write_seqs};

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

    // Assume all barcodes have same length
    let Some(&first_key) = barcode_data.keys().next() else {
        return Err(anyhow::anyhow!("Barcode data is empty"));
    };
    let bc_len = first_key.len();

    // Cache barcode keys (avoid repeated hashmap lookups)
    let barcodes: Vec<&[u8]> = barcode_data.keys().copied().collect();

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

        let matched = barcodes
            .iter()
            .find(|&&bc| bc_cmp(bc, &record.seq()[..bc_len], mismatch));

        match matched {
            Some(&bc) => {
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

    let (reverse_reader, _compression) = niffler::send::from_path(reverse)?;

    // Get records
    let mut forward_fastx_reader = needletail::parse_fastx_reader(forward_reader)?;
    //forward_records = forward_records.records();
    let mut reverse_fastx_reader = needletail::parse_fastx_reader(reverse_reader)?;

    // Clone barcode values in barcode_data structure for future iteration
    let my_vec = barcode_data.keys().cloned().collect::<Vec<_>>();

    // Get barcode length
    let bc_len = my_vec[0].len();

    // Initialize unknown files as empty
    let mut unk1_empty = "true";
    let mut unk2_empty = "true";

    // Change output compression format to user wanted compression
    // format if specified by --format option
    if format != niffler::send::compression::Format::No {
        compression = format;
    }

    while let Some(r) = forward_fastx_reader.next() {
        let record = r.expect("invalid record");
        let mut iter = my_vec.iter();
        let matched_barcode = iter.find(|&&x| bc_cmp(x, &record.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_seqs(
                &barcode_data.get(i).unwrap()[0],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        } else {
            unk1_empty = "false";
            write_seqs(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        }
    }

    while let Some(r) = reverse_fastx_reader.next() {
        let record = r.expect("invalid record");
        let mut iter = my_vec.iter();
        let matched_barcode = iter.find(|&&x| bc_cmp(x, &record.seq()[..bc_len], mismatch));

        if let Some(i) = matched_barcode {
            nb_records.entry(i).and_modify(|e| *e += 1).or_insert(1);
            write_seqs(
                &barcode_data.get(i).unwrap()[1],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        } else {
            unk2_empty = "false";
            write_seqs(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[1],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        }
    }
    let mut final_str = String::with_capacity(unk1_empty.len() + unk2_empty.len());
    final_str.push_str(unk1_empty);
    final_str.push_str(unk2_empty);

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
