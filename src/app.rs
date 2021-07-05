use clap::{crate_version, App, AppSettings, Arg};

pub fn build_app() -> App<'static, 'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        AppSettings::ColoredHelp
    } else {
        AppSettings::ColorNever
    };

    let app = App::new("sabreur")
        .version(crate_version!())
        .usage("sabreur [FLAGS/OPTIONS] <BARCODE> <FORWARD FILE> [<REVERSE FILE>]")
        .setting(clap_color_setting)
        .setting(AppSettings::DeriveDisplayOrder)
        .after_help(
            "Note: `sabreur -h` prints a short and concise overview while `sabreur --help` gives all \
                 details.",
        )
        .author("Anicet Ebou, anicet.ebou@gmail.com")
        .about("Fast, reliable and handy barcode demultiplexing tool for fastx files")
        .arg(
            Arg::with_name("BARCODE")
                .help("Input barcode file.")
                .long_help("Takes the barcode file containing barcode and output files data.\n \
                        Barcode file is tsv formated:\n \
                         `barcode1  file2_R1.fq  file1_R2.fq`\n \
                         `barcode2  file2_R1.fq  file2_R2.fq`\n \
                         `...`\n \
                        for paired-end data or like:\n \
                         `barcode1  file1.fq`\n \
                         `barcode2  file2.fq`\n \
                         `...`\n \
                        for single-end data.",
                )
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("FORWARD")
                .help("Input forward fasta or fastq file.\n Can be .(gz|xz|bz2) compressed.")
                .long_help(
                    "Input fasta or fastq forward file if demultiplexing paired-end\n \
                        data or to the single file in demultiplexing single-end data.",
                )
                .required(true)
                .index(2),
        )
        .arg(
            Arg::with_name("REVERSE")
                .help("Input reverse fasta or fastq file.\n Can be .(gz|xz|bz2) compressed.")
                .long_help(
                    "Input fasta or fastq reverse file if demultiplexing paired-end\n \
                        data. Should be ommited in single-end mode.",
                )
                .index(3),
        )
        .arg(
            Arg::with_name("mismatch")
                .long_help("Maximum number of mismatches allowed in barcode")
                .short("m")
                .long("mismatch")
                .value_name("N")
                .default_value("0"),
        )
        .arg(
            Arg::with_name("output")
                .help("Specifies the ouput directory")
                .short("o")
                .long("out")
                .value_name("FOLDER")
                .default_value("sabreur_out"),
        )
        .arg(
            Arg::with_name("format")
                .help("Set output files compression format.")
                .long_help(
                    "Specifies the compression format of the demultiplexed files:\n \
                        gz: for gzip files\n \
                        xz: for xz (lzma) files\n \
                        bz2: for bzip2 files\n \
                    Note: These options are available depending on your\n \
                          installation of their supporting libraries.\n \
                          Find more on sabreur homepage.",
                )
                .long("format")
                .short("f")
                .takes_value(true)
                .possible_values(&["gz", "xz", "bz2"])
                .hide_possible_values(true),
        )
        .arg(
            Arg::with_name("level")
                .help("Set the compression level")
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
                .short("l")
                .takes_value(true)
                .possible_values(&["1", "2", "3", "4", "5", "6", "7", "8", "9"])
                .hide_possible_values(true)
                .default_value("1"),
        )
        .arg(
            Arg::with_name("force")
                .help("Force reuse of default output directory")
                .long_help(
                    "Reuse the default output directory (sabreur_out).\n \
                    This will erase existing directory before creating it.",
                )
                .long("force")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("quiet")
                .long_help("Decrease program verbosity")
                .short("q")
                .long("quiet")
                .takes_value(false),
        );

    app
}
