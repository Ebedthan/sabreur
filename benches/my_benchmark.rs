use anyhow::Result;
use bio::io::fasta;
use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::Throughput;
use criterion::{black_box, criterion_group, criterion_main};
use std::fs::OpenOptions;
use std::io::{self};

fn bc_cmp(bc: &[u8]) -> bool {
    //let seq = "ATCGGCATATGCGATGCAGTAGTCAGTATGCAAAAACCCGGTTGGTTCCAA";
    //let s = &seq[..bc.len()];
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

fn write_to_fa<'a>(buffer: usize) -> Result<()> {
    let record =
        fasta::Record::with_attrs("id_str", Some("desc"), b"ATCGCCG");
    let compression = niffler::compression::Format::Gzip;
    let filename = "./myfile.fa";
    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(filename)?;

    let raw_out = Box::new(io::BufWriter::with_capacity(buffer, file));

    let handle = niffler::get_writer(
        raw_out,
        compression,
        niffler::compression::Level::One,
    )?;

    let mut writer = fasta::Writer::new(handle);
    let _fin_writer = writer.write_record(&record)?;

    Ok(())
}

fn from_write(c: &mut Criterion) {
    static KB: usize = 1024;

    let mut group = c.benchmark_group("from_write");
    for size in [KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            size,
            |b, &size| {
                b.iter(|| write_to_fa(size));
            },
        );
    }
    group.finish();
}

criterion_group!(benches, from_elem, from_write);
criterion_main!(benches);
