use criterion::Criterion;
use criterion::{black_box, criterion_group, criterion_main};
use triple_accel::*;

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
    c.bench_function("bc cmp", |b| {
        b.iter(|| bc_cmp(black_box(b"ATCGATGTGC"), black_box(b"ATCGATGTGA")))
    });
}

fn triple_bench(c: &mut Criterion) {
    c.bench_function("triple", |b| {
        b.iter(|| bc_triple(black_box(b"ATCGATGTGC"), black_box(b"ATCGATGTGA")))
    });
}

criterion_group!(benches, cmp_bench, triple_bench,);

criterion_main!(benches);
