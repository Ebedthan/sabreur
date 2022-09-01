// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;

use anyhow::Result;
use bio::io::{fasta, fastq};

use crate::io::{read_file, write_fa, write_fq};
use crate::utils::bc_cmp;

pub type Barcode<'a> = HashMap<&'a [u8], Vec<std::fs::File>>;

// Demultiplex a fasta::Record of single-end file
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
            .expect("file name should be available");
        } else {
            is_unk_empty = false;
            write_fa(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        }
    }
    Ok((nb_records, is_unk_empty))
}

// Demultiplex a fastq::Record of single-end file
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
            .expect("file name should be available");
        } else {
            is_unk_empty = false;
            write_fq(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &record,
                level,
            )
            .expect("file name should be available");
        }
    }
    Ok((nb_records, is_unk_empty))
}

// Demultiplex a fasta::Record of paired-end file
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
            .expect("file name should be available");
        } else {
            unk1_empty = "false";
            write_fa(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &f_rec,
                level,
            )
            .expect("file name should be available");
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
            .expect("file name should be available");
        } else {
            unk2_empty = "false";
            write_fa(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[1],
                compression,
                &r_rec,
                level,
            )
            .expect("file name should be available");
        }
    }
    let final_str = format!("{}{}", unk1_empty, unk2_empty);
    Ok((nb_records, final_str))
}

// Demultiplex a fasta::Record of paired-end file
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
            .expect("file name should be available");
        } else {
            unk1_empty = "false";
            write_fq(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[0],
                compression,
                &f_rec,
                level,
            )
            .expect("file name should be available");
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
            .expect("file name should be available");
        } else {
            unk2_empty = "false";
            write_fq(
                &barcode_data.get(&"XXX".as_bytes()).unwrap()[1],
                compression,
                &r_rec,
                level,
            )
            .expect("file name should be available");
        }
    }
    let final_str = format!("{}{}", unk1_empty, unk2_empty);
    Ok((nb_records, final_str))
}

// Tests ----------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

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
}
