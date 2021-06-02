use anyhow::Result;
use bio::io::fasta;
use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::{black_box, criterion_group, criterion_main};

fn bc_cmp(bc: &[u8]) -> bool {
    let mismatch = 0;
    let s = b"ATCG";

    // This wonderful line below compute the number of
    // character mismatch between two strings
    let nb_mismatch: i32 =
        bc.iter().zip(s.iter()).map(|(a, b)| (a != b) as i32).sum();

    nb_mismatch <= mismatch
}

fn from_elem(c: &mut Criterion) {
    let bc = b"ATCG";

    c.bench_with_input(
        BenchmarkId::new("bc_cmp_example", "ATCG"),
        &bc,
        |b, &x| {
            b.iter(|| bc_cmp(black_box(x)));
        },
    );
}

fn split_seq(record: &fasta::Record) {
    &record.seq()[..4];
}

fn from_split_seq(c: &mut Criterion) {
    let rec = fasta::Record::with_attrs("id1", Some("desc"), b"ATCGATCATCGATCGATCGAGTTGGTTTTGGGGTTTGGGTTTCCCAAAATGTTGATGTGTTTGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTT");

    c.bench_with_input(BenchmarkId::new("split_seq", "ATCG"), &rec, |b, x| {
        b.iter(|| split_seq(black_box(x)));
    });
}

fn write_to_fa<'a>(file: &'a std::fs::File) -> Result<()> {
    let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
    let cmp = niffler::compression::Format::Gzip;
    let handle = niffler::get_writer(
        Box::new(file),
        cmp,
        niffler::compression::Level::One,
    )?;

    let mut writer = fasta::Writer::new(handle);
    let _write_res = writer.write_record(&record)?;

    Ok(())
}

fn from_write(c: &mut Criterion) {
    let file = tempfile::tempfile().expect("Cannot create temp file");

    c.bench_with_input(
        BenchmarkId::new("from_write", "fasta"),
        &file,
        |b, file| {
            b.iter(|| write_to_fa(black_box(file)));
        },
    );
}

fn write_to_fa_sp<'a>(file: &'a std::fs::File) -> Result<()> {
    let cmp = niffler::compression::Format::Gzip;
    let handle = niffler::get_writer(
        Box::new(file),
        cmp,
        niffler::compression::Level::One,
    )?;

    let n_rec = fasta::Record::with_attrs("id", Some("desc"), &b"ATCGATCATCGATCGATCGAGTTGGTTTTGGGGTTTGGGTTTCCCAAAATGTTGATGTGTTTGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTT"[4..]);

    let mut writer = fasta::Writer::new(handle);
    let _write_res = writer.write_record(&n_rec)?;

    Ok(())
}

fn from_write_sp(c: &mut Criterion) {
    let file = tempfile::tempfile().expect("Cannot create temp file");

    c.bench_with_input(
        BenchmarkId::new("from_write_sp", "fasta"),
        &file,
        |b, file| {
            b.iter(|| write_to_fa_sp(black_box(file)));
        },
    );
}

criterion_group!(
    benches,
    from_elem,
    from_split_seq,
    from_write,
    from_write_sp
);
criterion_main!(benches);
