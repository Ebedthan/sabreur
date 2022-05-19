// Copyright 2021-2022 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use crate::utils;
use thiserror::Error;

use std::string::String;

#[derive(Debug, Error)]
pub enum Error {
    #[error(
        "Could not read file '{filename:}'. Does it exist and can be read by the user?"
    )]
    CantReadFile { filename: String },

    #[error("Could not create/read file '{filename:}'. Does the directory in path exist or can be written by the user?")]
    CantWriteFile { filename: String },

    #[error("Could not detect file format of '{filename:}'. File name should contains .fasta, .fa, .fastq, or .fq.")]
    UnableToDetectFileFormat { filename: String },

    #[error(
        "Error during writing of sample/species file in format {format:?}"
    )]
    WritingErrorNoFilename { format: utils::FileType },
}
