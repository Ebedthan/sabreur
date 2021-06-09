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

fn to_niffler_level_match(int_level: i32) -> niffler::Level {
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

fn to_niffler_level_if(int_level: i32) -> niffler::Level {
    if int_level == 1 {
        niffler::Level::One
    } else if int_level == 2 {
        niffler::Level::Two
    } else if int_level == 3 {
        niffler::Level::Three
    } else if int_level == 4 {
        niffler::Level::Four
    } else if int_level == 5 {
        niffler::Level::Five
    } else if int_level == 6 {
        niffler::Level::Six
    } else if int_level == 7 {
        niffler::Level::Seven
    } else if int_level == 8 {
        niffler::Level::Eight
    } else if int_level == 9 {
        niffler::Level::Nine
    } else {
        niffler::Level::One
    }
}

fn from_level_match(c: &mut Criterion) {
    c.bench_function("from_level_match", |b| {
        b.iter(|| to_niffler_level_match(black_box(2)))
    });
}

fn from_level_if(c: &mut Criterion) {
    c.bench_function("from_level_match", |b| {
        b.iter(|| to_niffler_level_if(black_box(2)))
    });
}

criterion_group!(
    benches,
    from_elem,
    from_split_seq,
    from_write,
    from_level_match,
    from_level_if,
);
criterion_main!(benches);
