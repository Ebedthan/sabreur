# sabreur
[![Continuous Integration](https://github.com/Ebedthan/sabreur/actions/workflows/ci.yml/badge.svg)](https://github.com/Ebedthan/sabreur/actions/workflows/ci.yml)
<a href="https://crates.io/crates/sabreur">
    <img src="https://img.shields.io/crates/v/sabreur.svg?style=flat">
</a>
<a href="https://codecov.io/gh/Ebedthan/sabreur">
    <img src="https://codecov.io/gh/Ebedthan/sabreur/branch/main/graph/badge.svg">
</a>
<a href="https://github.com/Ebedthan/sabreur">
    <img src="https://img.shields.io/badge/rust-1.74.1%2B-blue.svg?maxAge=3600">
</a>
<a href="https://github.com/Ebedthan/sabreur/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue?style=flat">
</a>

# About

With next-generation sequencing tools capabilities, millions to billions of reads are generated. To reach such a rate in a cost-efficient manner, barcoding individual sequences for multiple lines or species is a common practice.

Sabreur is a tool that aims to demultiplex barcoded reads into separate files. It supports both fasta and fastq files. Input files can be gzip, bzip2 or xz compressed in input or output (Thanks to the awesome [niffler crate](https://github.com/luizirber/niffler)). If an uncompressed file is provided the output is by default uncompressed. But this behaviour can be changed by settingn the `--format` option to the desired compress format. The `--format` option if specified while input files are compressed changes output files to the specified compress format. Sabreur in its core compares the provided barcodes with each read, then separates the read into its appropriate file. If a read does not have a recognized barcode, then it is put into an unknown file.


## How to use sabreur

### Paired-end mode
```
sabreur barcode.txt input_R1.fq.gz input_R2.fq.gz
```

### Single-end mode
```
sabreur barcode.txt input.fq
```

Input sequences files can be fasta or fastq, compressed or not. 
The supported compression format are gz, bz2, xz and zst.
Just give the sequences, sabreur know how to handle it!

## Command-line arguments

```
USAGE:
    sabreur [options] <BARCODE> <FORWARD FILE> [<REVERSE FILE>]

ARGS:
    <BARCODE>    input barcode file
    <FORWARD>    input forward fastx file
    <REVERSE>    input reverse fastx file

OPTIONS:
    -m, --mismatch <INT>    maximum number of mismatches [default: 0]
    -o, --out <DIR>         ouput directory [default: sabreur_out]
    -f, --format <STR>      output files compression format
    -l, --level <INT>       compression level [default: 1]
        --force             force reuse of output directory
    -q, --quiet             decrease program verbosity
    -h, --help              Print help information
    -V, --version           Print version information

```

## Requirements
- [Rust](https://rust-lang.org) in stable channel
- libgz for gz file support
- liblzma for xz file support
- libbzip2 for bzip2 file support
- zstd for zstd file support


## Installation

## From crates.io
If you already have a functional rust installation do:

```
cargo install sabreur
```

## From source
```
git clone https://github.com/Ebedthan/sabreur.git
cd sabreur

cargo build --release
cargo test
cargo install --path .
```

## Benchmark

We used [hyperfine](https://github.com/sharkdp/hyperfine) for benchmarking with this [dataset](https://figshare.com/articles/dataset/Paired-end_fastq_files_for_demultiplexing/14701629).


| Tool  | Single-end uncompressed output | Single-end compressed output | Paired-end uncompressed output | Paired-end compressed output |
| :---  |             :----:             |             :----:           |              :----:           |              :----:           |
| [idemp](https://github.com/yhwu/idemp) | - | 211.571 ± 3.718 | -      | 366.247 ± 10.482  |
| [sabre](https://github.com/najoshi/sabre) | 32.911 ± 2.411 | - | 109.470 ± 49.909 | -     |
| **sabreur** | 10.843 ± 0.531| 93.840 ± 0.446    | 40.878 ± 13.743     | 187.533 ± 0.572   |



A simple benchmark of the different compression format, zst being the fastest.
| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `sabreur tests/bc_pe_fq.txt tests/input_R1.fastq.gz tests/input_R2.fastq.gz --format zst` | 43.096 ± 1.547 | 41.179 | 46.878 | 1.00 |
| `sabreur tests/bc_pe_fq.txt tests/input_R1.fastq.gz tests/input_R2.fastq.gz --format bz2` | 94.049 ± 4.762 | 87.984 | 101.140 | 2.18 ± 0.14 |
| `sabreur tests/bc_pe_fq.txt tests/input_R1.fastq.gz tests/input_R2.fastq.gz (--format gz)` | 123.107 ± 1.748 | 120.529 | 125.166 | 2.86 ± 0.11 |
| `sabreur tests/bc_pe_fq.txt tests/input_R1.fastq.gz tests/input_R2.fastq.gz --format xz` | 285.692 ± 18.625 | 264.960 | 325.750 | 6.63 ± 0.49 |


## Note
Sabreur use colored output in help, nevertheless sabreur honors [NO_COLORS](https://no-color.org/) environment variable.

Sabreur use a special barcode tab-delimited file format in the form:

```
barcode1    barcode1_file1.fq   barcode1_file2.fq
barcode2    barcode2_file1.fq   barcode2_file2.fq
...
```

### Minimum supported Rust version
`sabreur` minimum [Rust](https://www.rust-lang.org/) version is 1.74.1.

## Contributions
Contributions are welcomed under the project [code of conduct](https://github.com/Ebedthan/sabreur#code-of-conduct).

## Bugs
Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/sabreur/issues).

## Code of conduct
Please note that the sabreur project is released with a [Contributor Code of Conduct](https://github.com/Ebedthan/sabreur/blob/main/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.


<p align="center">
    <a href="https://github.com/Ebedthan/sabreur">
        <img src="img/sabreur.png" width="300">
    </a>
</p>
