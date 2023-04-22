// Copyright 2021-2023 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{crate_version, AppSettings, Arg, ColorChoice, Command};

pub fn build_app() -> Command<'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        ColorChoice::Always
    } else {
        ColorChoice::Never
    };

    let app = Command::new("sabreur")
        .version(crate_version!())
        .override_usage("sabreur [options] <BARCODE> <FORWARD FILE> [<REVERSE FILE>]")
        .color(clap_color_setting)
        .setting(AppSettings::DeriveDisplayOrder)
        .after_help(
            "Note: `sabreur -h` prints a short and concise overview while `sabreur --help` gives all \
                 details.",
        )
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Fast, reliable and handy barcode demultiplexing for fastx files")
        .arg(
            Arg::new("BARCODE")
                .help("input barcode file")
                .long_help("Takes the barcode file containing barcode and output files data\n \
                        Barcode file is tsv formated:\n \
                         `barcode1  file2_R1.fq  file1_R2.fq`\n \
                         `barcode2  file2_R1.fq  file2_R2.fq`\n \
                         `...`\n \
                        for paired-end data or like:\n \
                         `barcode1  file1.fq`\n \
                         `barcode2  file2.fq`\n \
                         `...`\n \
                        for single-end data",
                )
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("FORWARD")
                .help("input forward fastx file\n")
                .long_help(
                    "Input fasta or fastq forward file if demultiplexing paired-end\n \
                        data or to the single file in demultiplexing single-end data",
                )
                .required(true)
                .index(2),
        )
        .arg(
            Arg::new("REVERSE")
                .help("input reverse fastx file\n")
                .long_help(
                    "Input fasta or fastq reverse file if demultiplexing paired-end\n \
                        data. Should be ommited in single-end mode",
                )
                .index(3),
        )
        .arg(
            Arg::new("mismatch")
                .help("maximum number of mismatches")
                .long_help("maximum number of mismatches allowed in a barcode ")
                .short('m')
                .long("mismatch")
                .value_name("INT")
                .default_value("0"),
        )
        .arg(
            Arg::new("output")
                .help("ouput directory")
                .short('o')
                .long("out")
                .value_name("DIR")
                .default_value("sabreur_out"),
        )
        .arg(
            Arg::new("format")
                .help("output files compression format")
                .long_help(
                    "Specifies the compression format of the demultiplexed files:\n \
                        gz: for gzip files\n \
                        xz: for xz (lzma) files\n \
                        bz2: for bzip2 files\n \
                        zst: for zstd files \n \
                    Note: These options are available depending on your\n \
                          installation of their supporting libraries.\n \
                          Find more on sabreur homepage",
                )
                .long("format")
                .short('f')
                .takes_value(true)
                .value_name("STR")
                .possible_values(["gz", "xz", "bz2", "zst"])
                .hide_possible_values(true),
        )
        .arg(
            Arg::new("level")
                .help("compression level")
                .long_help(
                    "Specifies the compression level wanted for the demultiplexed file:\n \
                        1: Level One, optimize the compression time\n \
                        2: Level Two\n \
                        3: Level Three\n \
                        4: Level Four\n \
                        5: Level Five\n \
                        6: Level Six\n \
                        7: Level Seven\n \
                        8: Level Eight\n \
                        9: Level Nine, optimize the size of the output\n",
                )
                .long("level")
                .short('l')
                .takes_value(true)
                .value_name("INT")
                .possible_values(["1", "2", "3", "4", "5", "6", "7", "8", "9"])
                .hide_possible_values(true)
                .default_value("1"),
        )
        .arg(
            Arg::new("force")
                .help("force reuse of output directory")
                .long_help(
                    "Reuse the default output directory (sabreur_out).\n \
                    This will erase existing directory before creating it.",
                )
                .long("force")
                .takes_value(false),
        )
        .arg(
            Arg::new("quiet")
                .long_help("decrease program verbosity")
                .short('q')
                .long("quiet")
                .takes_value(false),
        );

    app
}
