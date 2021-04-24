[![License](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](https://github.com/Ebedthan/sabreur/blob/master/LICENSE)
![CI](https://github.com/Ebedthan/sabreur/workflows/CI/badge.svg)
[![CodeCov](https://codecov.io/gh/Ebedthan/sabreur/branch/main/graph/badge.svg)](https://codecov.io/gh/Ebedthan/sabreur)

# <img src="./img/ninja.png" width=40em alt="sabreur" /> 🧬 sabreur, a barcode demultiplexing tool for fasta and fastq files.

With next-generation sequencing tools capabilities, millions to billions of reads are generated. To reach such a rate in a cost-efficient manner, barcoding individual sequences for multiple lines or species is a common practice.

Sabreur is a tool that aims to demultiplex barcoded reads into separate files. It supports both fasta and fastq files either gzipped or not. The resulting files are compressed or not following the compression mode of input files. Sabreur in its core compares the provided barcodes with each read, then separates the read into its appropriate file. If a read does not have a recognized barcode, then it is put into an unknown file.


## How to use sabreur

### Paired-end mode
```
sabreur barcode.txt input_R1.fq.gz input_R2.fq.gz
```

### Single-end mode
```
sabreur barcode.txt input.fa.gz
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
    -m, --mismatch <N>    Maximum number of mismatches allowed in barcode [default: 0]
    -o, --out <FOLDER>    Output folder [default: sabreur_out]

ARGS:
    <BARCODE>    Input barcode file [required]
    <FORWARD>    Input forward fasta or fastq file. Can be gzipped [required]
    <REVERSE>    Input reverse fasta or fastq file. Can be gzipped
```

## Requirements
- [Rust](https://rust-lang.org) in stable channel
- libgz


## Installation

### From source
```
git clone https://github.com/Ebedthan/sabreur.git
cd sabreur

cargo build
cargo test
cargo install --path .
```

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