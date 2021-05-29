[![Crates.io](https://img.shields.io/crates/v/sabreur.svg)](https://crates.io/crates/sabreur)
[![License](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](https://github.com/Ebedthan/sabreur/blob/master/LICENSE)
![CI](https://github.com/Ebedthan/sabreur/workflows/CI/badge.svg)
[![CodeCov](https://codecov.io/gh/Ebedthan/sabreur/branch/main/graph/badge.svg)](https://codecov.io/gh/Ebedthan/sabreur)

# sabreur: fast, reliable and handy demultiplexing tool for fastx files.

With next-generation sequencing tools capabilities, millions to billions of reads are generated. To reach such a rate in a cost-efficient manner, barcoding individual sequences for multiple lines or species is a common practice.

Sabreur is a tool that aims to demultiplex barcoded reads into separate files. It supports both fasta and fastq files. Input files can be gzip, bzip2 or xz compressed in input or output (Thanks to the awesome [niffler crate](https://github.com/luizirber/niffler)). If an uncompressed file is provided the output is by default uncompressed. But this behaviour can be changed by settingn the `--format` option to the desired compress format. The `--format` option if specified while input files are compressed changes output files to the specified compress format. Sabreur in its core compares the provided barcodes with each read, then separates the read into its appropriate file. If a read does not have a recognized barcode, then it is put into an unknown file.


## How to use sabreur

### Paired-end mode
```
sabreur barcode.txt input_R1.fq.gz input_R2.fq.gz
```

### Single-end mode
```
sabreur barcode.txt input.fa --format xz
```

Input sequences files can be fasta or fastq, gzipped or not. Just give the sequences, sabreur know how to handle it!

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

 

## Note
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

## License
Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.

## Code of conduct
Please note that the sabreur project is released with a [Contributor Code of Conduct](https://github.com/Ebedthan/sabreur/blob/main/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.