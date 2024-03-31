// Copyright 2021-2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs::File;
use std::io;

use anyhow::Result;

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
) -> Result<()> {
    let mut handle = niffler::send::get_writer(Box::new(file), compression, level)?;

    match record.format() {
        needletail::parser::Format::Fasta => needletail::parser::write_fasta(record.id(), &record.seq(), &mut handle, needletail::parser::LineEnding::Unix)?,
        needletail::parser::Format::Fastq => needletail::parser::write_fastq(record.id(), &record.seq(), record.qual(), &mut handle, needletail::parser::LineEnding::Unix)?
    }
    
    Ok(())
}

// Tests --------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

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
