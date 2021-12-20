use anyhow::Result;
use bio::io::fasta;
use criterion::Criterion;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Throughput};
use std::fs::{self, File};
use std::io;
use triple_accel::*;

fn write_fa_capacity<'a>(
    file: &'a std::fs::File,
    capacity: usize,
) -> Result<()> {
    let record = fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
    let cmp = niffler::compression::Format::Gzip;
    let handle = niffler::get_writer(
        Box::new(file),
        cmp,
        niffler::compression::Level::One,
    )?;

    let mut writer = fasta::Writer::with_capacity(capacity, handle);
    let _write_res = writer.write_record(&record)?;

    Ok(())
}

fn write_fa_capacity_benchmarck(c: &mut Criterion) {
    let myfile = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open("tests/write_fa_default.fa")
        .expect("Cannot open file");

    static KB: usize = 8192;

    let mut group = c.benchmark_group("write_file");
    for size in [KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            size,
            |b, &size| {
                b.iter(|| write_fa_capacity(&myfile, size));
            },
        );
    }
    group.finish();
}

fn read_file_capacity(filename: &str, capacity: usize) -> Result<()> {
    let raw_in = Box::new(io::BufReader::new(
        File::open(filename).expect("Cannot open file"),
    ));

    let (reader, _compression) =
        niffler::get_reader(raw_in).expect("Cannot open file");

    let _records = fasta::Reader::with_capacity(capacity, reader).records();

    Ok(())
}

fn read_file_capacity_benchmarck(c: &mut Criterion) {
    static KB: usize = 8192;

    let mut group = c.benchmark_group("read_file");
    for size in [KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            size,
            |b, &size| {
                b.iter(|| read_file_capacity("tests/input_R1.fastq.gz", size));
            },
        );
    }
    group.finish();
}

fn bc_cmp(bc: &[u8], seq: &[u8]) -> bool {
    // This wonderful line below compute the number of
    // character mismatch between two strings
    let mismatch = 1;
    let nb_mismatch: i32 = bc
        .iter()
        .zip(seq.iter())
        .map(|(a, b)| (a != b) as i32)
        .sum();

    nb_mismatch <= mismatch
}

fn bc_triple(bc: &[u8], seq: &[u8]) -> bool {
    hamming(bc, seq) <= 1
}

fn cmp_bench(c: &mut Criterion) {
    c.bench_function("bc cmp", |b| b.iter(|| bc_cmp(black_box(b"ATCGATGTGC"), black_box(b"ATCGATGTGA"))));
}

fn triple_bench(c: &mut Criterion) {
    c.bench_function("triple", |b| b.iter(|| bc_triple(black_box(b"ATCGATGTGC"), black_box(b"ATCGATGTGA"))));
}

criterion_group!(
    benches,
    cmp_bench,
    triple_bench,
);

criterion_main!(benches);
