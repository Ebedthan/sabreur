<p align="center">
    <a href="https://github.com/Ebedthan/sabreur">
        <img src="img/sabreur.png" width="300">
    </a>
    </br>
    <a href="https://github.com/Ebedthan/sabreur/actions?query=workflow%3A%22Continuous+Integration%22">
        <img src="https://img.shields.io/github/workflow/status/Ebedthan/sabreur/Continuous%20Integration?style=flat&logo=GitHub%20Actions">
    </a>
    <a href="https://github.com/Ebedthan/sabreur/actions?query=workflow%3A%22Continuous+Deployment%22">
        <img src="https://img.shields.io/github/workflow/status/Ebedthan/sabreur/Continuous%20Deployment?style=flat&logo=GitHub%20Actions&label=deploy">
    </a>
    <a href="https://crates.io/crates/sabreur">
        <img src="https://img.shields.io/crates/v/sabreur.svg?style=flat">
    </a>
    <a href="https://codecov.io/gh/Ebedthan/sabreur">
        <img src="https://codecov.io/gh/Ebedthan/sabreur/branch/main/graph/badge.svg">
    </a>
    [![Rust](https://img.shields.io/badge/rust-1.56.1%2B-blue.svg?maxAge=3600)](https://github.com/Ebedthan/sabreur)
    <a href="https://github.com/Ebedthan/sabreur/blob/master/LICENSE">
        <img src="https://img.shields.io/badge/license-MIT-blue?style=flat">
    </a>
</p>

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
    sabreur [FLAGS] [OPTIONS] <BARCODE> <FORWARD> [REVERSE]

FLAGS:
        --force      Force reuse of output directory
    -h, --help       Prints help information
    -q, --quiet      Decrease program verbosity
    -V, --version    Prints version information

OPTIONS:
    -f, --format <format>    Set output files compression format. [possible values: gz, xz, bz2]
    -l, --level <level>      Set the compression level [default: 1]  [possible values: 1, 2, 3, 4, 5, 6, 7, 8, 9]
    -m, --mismatch <N>       Maximum number of mismatches allowed in a barcode [default: 0]
    -o, --out <FOLDER>       Output folder [default: sabreur_out]

ARGS:
    <BARCODE>    Input barcode file.
    <FORWARD>    Input forward fasta or fastq file. Can be gz, xz or bz2 compressed.
    <REVERSE>    Input reverse fasta or fastq file. Can be gz, xz or bz2 compressed.
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


### Minimum Rust version policy
This crate's minimum supported `rustc` version is `1.56.1`.

## Note
Sabreur use colored output in help, nevertheless sabreur honors [NO_COLORS](https://no-color.org/) environment variable.

Sabreur use a special barcode tab-delimited file format in the form:

```
barcode1    barcode1_file1.fq   barcode1_file2.fq
barcode2    barcode2_file1.fq   barcode2_file2.fq
...
```

## Contributions
Contributions are welcomed under the project [code of conduct](https://github.com/Ebedthan/sabreur#code-of-conduct).

## Bugs
Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/sabreur/issues).

## Code of conduct
Please note that the sabreur project is released with a [Contributor Code of Conduct](https://github.com/Ebedthan/sabreur/blob/main/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.