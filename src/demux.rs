// Copyright 2021-2026 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;
use std::io::Write;

use crate::utils::{self, write_seqs};

pub type Barcode<'a> = HashMap<&'a [u8], Vec<Box<dyn Write + Send + 'static>>>;

/// A function to demultiplex a FASTA/FASTQ file
pub fn se_demux<'a, 'n>(
    file: &str,
    barcode_data: &mut Barcode<'a>,
    mismatch: u8,
    nb_records: &'n mut HashMap<&'a [u8], u32>,
) -> anyhow::Result<(&'n mut HashMap<&'a [u8], u32>, bool)> {
    // Prepare decompression stream (input decompression is auto-detected
    // from the file itself; output compression was already decided when
    // each writer in `barcode_data` was created, so it isn't needed here)
    let (reader, _) = niffler::send::from_path(file)?;
    let mut reader = needletail::parse_fastx_reader(reader)?;

    // Cache barcode keys (avoid repeated hashmap lookups). The "XXX"
    // sentinel is excluded: it is a fallback output slot, not a real
    // barcode, and must never be a match candidate. Sorting makes
    // ambiguous-match tie-breaking (under mismatch > 0) deterministic
    // across runs, instead of depending on HashMap iteration order.
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

    // Process each record
    while let Some(record) = reader.next() {
        let record = record?;

        let seq = record.seq();
        let matched = (seq.len() >= bc_len)
            .then(|| utils::best_barcode_match(&barcodes, &seq[..bc_len], mismatch))
            .flatten();

        match matched {
            Some(bc) => {
                *nb_records.entry(bc).or_insert(0) += 1;
                let files = barcode_data
                    .get_mut(bc)
                    .expect("matched barcode must exist in barcode_data");
                write_seqs(files[0].as_mut(), &record)?;
            }
            None => {
                is_unk_empty = false;
                let files = barcode_data.get_mut(b"XXX".as_ref()).ok_or_else(|| {
                    anyhow::anyhow!("Missing 'XXX' fallback barcode entry in barcode_data")
                })?;
                write_seqs(files[0].as_mut(), &record)?;
            }
        }
    }
    Ok((nb_records, is_unk_empty))
}

/// A function to demultiplex a pair of FASTA/FASTQ files
pub fn pe_demux<'a, 'n>(
    forward: &str,
    reverse: &str,
    barcode_data: &mut Barcode<'a>,
    mismatch: u8,
    nb_records: &'n mut HashMap<&'a [u8], u32>,
) -> anyhow::Result<(&'n mut HashMap<&'a [u8], u32>, String)> {
    // Get fasta files readers (input decompression is auto-detected from
    // each file; output compression was already decided when each writer
    // in `barcode_data` was created)
    let (forward_reader, _) = niffler::send::from_path(forward)?;
    let (reverse_reader, _) = niffler::send::from_path(reverse)?;

    let mut forward_fastx_reader = needletail::parse_fastx_reader(forward_reader)?;
    let mut reverse_fastx_reader = needletail::parse_fastx_reader(reverse_reader)?;

    // Get barcode information once. Exclude the "XXX" fallback slot from
    // the candidate list (it's an output destination, not a barcode), and
    // sort for deterministic matching order under mismatch > 0.
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

    // Process forward and reverse reads together, in lockstep, so that mates
    // always stay paired. The barcode is looked up from the forward (R1)
    // read only, and the resulting assignment (a barcode, or "unknown") is
    // applied to BOTH mates. This is what keeps output R1/R2 files
    // synchronized record-for-record.
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
                     records; they are not properly paired."
                ));
            }
        };

        let fwd_seq = fwd_record.seq();
        let matched = (fwd_seq.len() >= bc_len)
            .then(|| utils::best_barcode_match(&barcodes, &fwd_seq[..bc_len], mismatch))
            .flatten();

        match matched {
            Some(bc) => {
                *nb_records.entry(bc).or_insert(0) += 1;
                let files = barcode_data
                    .get_mut(bc)
                    .expect("matched barcode must exist in barcode_data");
                write_seqs(files[0].as_mut(), &fwd_record)?;
                write_seqs(files[1].as_mut(), &rev_record)?;
            }
            None => {
                unk1_empty = false;
                unk2_empty = false;
                let files = barcode_data.get_mut(b"XXX".as_ref()).ok_or_else(|| {
                    anyhow::anyhow!("Missing 'XXX' fallback barcode entry in barcode_data")
                })?;
                write_seqs(files[0].as_mut(), &fwd_record)?;
                write_seqs(files[1].as_mut(), &rev_record)?;
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

    // Helper: a tempfile boxed as a Write + Send trait object, matching
    // what utils::create_output_writer produces in the real binary.
    fn boxed_tempfile() -> Box<dyn Write + Send + 'static> {
        Box::new(tempfile::tempfile().expect("Cannot create temp file"))
    }

    #[test]
    fn test_se_demux_1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fa.gz", &mut bc_data, 0, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_trim() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fa.gz", &mut bc_data, 0, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_m1() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"ATTGTT", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fa.gz", &mut bc_data, 1, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_m2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"ATTGTT", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fa.gz", &mut bc_data, 2, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_2() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"ATTGTT", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fq.gz", &mut bc_data, 0, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_m3() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"ATTGTT", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fq.gz", &mut bc_data, 1, &mut nb_records).is_ok());
    }

    #[test]
    fn test_se_demux_m4() {
        let mut bc_data: Barcode = HashMap::new();
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        bc_data.insert(b"ACCGTA", vec![boxed_tempfile()]);
        bc_data.insert(b"ATTGTT", vec![boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);

        assert!(se_demux("tests/test.fq.gz", &mut bc_data, 2, &mut nb_records).is_ok());
    }

    // A named tempfile whose path we can reopen and read back after the
    // writer boxed into `Barcode` has been flushed. Returns (path, boxed
    // writer, guard). The guard must be kept alive for the duration of the
    // test (it deletes the file on drop).
    fn named_writer() -> (
        std::path::PathBuf,
        Box<dyn Write + Send + 'static>,
        tempfile::NamedTempFile,
    ) {
        let guard = tempfile::NamedTempFile::new().expect("Cannot create named temp file");
        let path = guard.path().to_path_buf();
        let writer: Box<dyn Write + Send + 'static> =
            Box::new(guard.reopen().expect("Cannot reopen temp file"));
        (path, writer, guard)
    }

    fn flush_all(bc_data: &mut Barcode) {
        for files in bc_data.values_mut() {
            for w in files.iter_mut() {
                w.flush().expect("flush failed");
            }
        }
    }

    #[test]
    fn test_se_demux_routes_reads_to_correct_files() {
        use std::io::Read;

        let (fwd_path, fwd_writer, _fwd_guard) = named_writer();
        let (unk_path, unk_writer, _unk_guard) = named_writer();

        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert(b"ACCGTA", vec![fwd_writer]);
        bc_data.insert(b"XXX", vec![unk_writer]);

        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        // tests/test.fa.gz contains: read1 (matches ACCGTA exactly),
        // read2 (ATTGTT-prefixed, no match since only ACCGTA is
        // registered here), read3 (no barcode match), read4 (shorter than
        // the barcode, must not panic, must land in unknown).
        let (stats, unk_empty) =
            se_demux("tests/test.fa.gz", &mut bc_data, 0, &mut nb_records).expect("demux failed");

        assert_eq!(stats.get(b"ACCGTA".as_ref()), Some(&1));
        assert!(!unk_empty);

        flush_all(&mut bc_data);

        let mut fwd_content = String::new();
        std::fs::File::open(&fwd_path)
            .unwrap()
            .read_to_string(&mut fwd_content)
            .unwrap();
        assert!(fwd_content.contains("read1"));
        assert!(!fwd_content.contains("read2"));
        assert!(!fwd_content.contains("read3"));
        assert!(!fwd_content.contains("read4"));

        let mut unk_content = String::new();
        std::fs::File::open(&unk_path)
            .unwrap()
            .read_to_string(&mut unk_content)
            .unwrap();
        assert!(unk_content.contains("read2"));
        assert!(unk_content.contains("read3"));
        assert!(unk_content.contains("read4"));
    }

    #[test]
    fn test_se_demux_empty_barcode_data_errors() {
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert(b"XXX", vec![boxed_tempfile()]);
        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        // Only the XXX sentinel is present, no real barcodes, so this
        // must fail with a clear error, not panic on an empty Vec.
        let result = se_demux("tests/test.fa.gz", &mut bc_data, 0, &mut nb_records);
        assert!(result.is_err());
    }

    #[test]
    fn test_pe_demux_keeps_mates_paired_and_routes_from_r1_only() {
        use std::io::Read;

        let (fwd_path, fwd_writer1, _g1) = named_writer();
        let (rev_path, fwd_writer2, _g2) = named_writer();
        let (unk_fwd_path, unk_writer1, _g3) = named_writer();
        let (unk_rev_path, unk_writer2, _g4) = named_writer();

        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert(b"ACCGTA", vec![fwd_writer1, fwd_writer2]);
        bc_data.insert(b"XXX", vec![unk_writer1, unk_writer2]);

        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        // test_R1.fq.gz / test_R2.fq.gz: pair1's R1 carries barcode
        // ACCGTA; its R2 mate does NOT start with ACCGTA. Both records of
        // pair1 must end up in the ACCGTA files together, because the
        // barcode decision is made from R1 only. pair2 matches nothing on
        // R1, so both of its mates go to unknown together.
        let (stats, unk_status) = pe_demux(
            "tests/test_R1.fq.gz",
            "tests/test_R2.fq.gz",
            &mut bc_data,
            0,
            &mut nb_records,
        )
        .expect("pe_demux failed");

        assert_eq!(stats.get(b"ACCGTA".as_ref()), Some(&1));
        assert_eq!(unk_status, "falsefalse");

        flush_all(&mut bc_data);

        let mut fwd_content = String::new();
        std::fs::File::open(&fwd_path)
            .unwrap()
            .read_to_string(&mut fwd_content)
            .unwrap();
        assert!(fwd_content.contains("pair1"));
        assert!(!fwd_content.contains("pair2"));

        let mut rev_content = String::new();
        std::fs::File::open(&rev_path)
            .unwrap()
            .read_to_string(&mut rev_content)
            .unwrap();
        // pair1's R2 mate must be here too, even though its own sequence
        // doesn't start with ACCGTA, it follows R1's assignment.
        assert!(rev_content.contains("pair1"));
        assert!(!rev_content.contains("pair2"));

        let mut unk_fwd_content = String::new();
        std::fs::File::open(&unk_fwd_path)
            .unwrap()
            .read_to_string(&mut unk_fwd_content)
            .unwrap();
        assert!(unk_fwd_content.contains("pair2"));

        let mut unk_rev_content = String::new();
        std::fs::File::open(&unk_rev_path)
            .unwrap()
            .read_to_string(&mut unk_rev_content)
            .unwrap();
        assert!(unk_rev_content.contains("pair2"));
    }

    #[test]
    fn test_pe_demux_errors_on_unequal_record_counts() {
        let mut bc_data: Barcode = HashMap::new();
        bc_data.insert(b"ACCGTA", vec![boxed_tempfile(), boxed_tempfile()]);
        bc_data.insert(b"XXX", vec![boxed_tempfile(), boxed_tempfile()]);

        let mut nb_records: HashMap<&[u8], u32> = HashMap::new();

        // test_R1.fq.gz has 2 records, short_R2.fq.gz has only 1, this
        // must be a clear error, not silently truncate to the shorter file.
        let result = pe_demux(
            "tests/test_R1.fq.gz",
            "tests/short_R2.fq.gz",
            &mut bc_data,
            0,
            &mut nb_records,
        );
        assert!(result.is_err());
    }
}
