// Copyright 2021-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "sabreur",
    author = "Anicet Ebou, anicet.ebou@gmail.com",
    version,
    about,
    override_usage = "sabreur [options] <BARCODE> <FORWARD FILE> [<REVERSE FILE>]",
    after_help = "Note: `sabreur -h` prints a short and concise overview while `sabreur --help` gives all details."
)]
pub struct Cli {
    /// Input barcode file
    #[arg(value_name = "BARCODE", value_parser = is_file)]
    pub barcode: String,

    /// Input forward fastx file
    #[arg(value_name = "FORWARD", value_parser = is_file)]
    pub forward: String,

    /// Input reverse fastx file (optional)
    #[arg(value_name = "REVERSE", value_parser = is_file)]
    pub reverse: Option<String>,

    /// Maximum number of mismatches
    #[arg(short, long, default_value_t = 0)]
    pub mismatch: u8,

    /// Output directory
    #[arg(short, long, default_value = "sabreur_out")]
    pub output: PathBuf,

    /// Output files compression format
    #[arg(short, long, value_enum, hide_possible_values = true)]
    pub format: Option<CompressionFormat>,

    /// Compression level (1-9)
    #[arg(short, long, default_value_t = 1, hide_possible_values = true)]
    pub level: u8,

    /// Force reuse of output directory
    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub force: bool,

    /// Decrease program verbosity
    #[arg(short, long, action = clap::ArgAction::SetTrue)]
    pub quiet: bool,
}

#[derive(Debug, Copy, Clone, ValueEnum, PartialEq, Eq)]
pub enum CompressionFormat {
    Gz,
    Xz,
    Bz2,
    Zst,
}

fn is_file(s: &str) -> Result<String, String> {
    if std::path::Path::new(s).is_file() {
        Ok(s.to_string())
    } else {
        Err("path does not exist".to_string())
    }
}
