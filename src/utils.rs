/*
Copyright (c) 2021 Anicet Ebou <anicet.ebou@gmail.com>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use anyhow::Result;


use std::io;
use std::path::Path;
// use std::fs::File;

/* ***********************
 * setup logging function
 * ***********************/
pub fn setup_logging(quiet: bool) -> Result<(), fern::InitError> {
    let mut base_config = fern::Dispatch::new();

    base_config = match quiet {
        true => base_config.level(log::LevelFilter::Info)
                        .level_for("overly-verbose-target", log::LevelFilter::Warn),
        _ => base_config.level(log::LevelFilter::Trace),
    };

    let file_config = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.target(),
                record.level(),
                message
            ))
        })
        .chain(fern::log_file("maph.log")?);
    
        let stdout_config = fern::Dispatch::new()
        .format(|out, message, record| {
            // special format for debug messages coming from maph.
            if record.level() > log::LevelFilter::Info && record.target() == "maph" {
                out.finish(format_args!(
                    "---\nDEBUG: {}: {}\n---",
                    chrono::Local::now().format("%H:%M:%S"),
                    message
                ))
            } else {
                out.finish(format_args!(
                    "[{}][{}] {}",
                    chrono::Local::now().format("%H:%M:%S"),
                    record.level(),
                    message
                ))
            }
        })
        .chain(io::stdout());

    base_config
        .chain(file_config)
        .chain(stdout_config)
        .apply()?;

    Ok(())
}

/* ***********************
 * read_gz function


pub fn read_file(p: &Path) -> Result<(Box<dyn io::Read>, niffler::compression::Format), > {
    let raw_in = Box::new(io::BufReader::new(
        File::open(p)?
    ));

    let (reader, compression) = niffler::get_reader(raw_in)
        .expect("Cannot read input fasta file");
    Ok((reader, compression))
    
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_gz() {
        let (mut reader, compression) = read_file(Path::new("tests/tgz.fa.gz")).unwrap();
        let mut contents = String::new();
        reader.read_to_string(&mut contents).expect("Error during file reading");

        assert_eq!(compression, niffler::compression::Format::Gzip);
        assert_eq!(contents, ">seqID1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n");
    }

    #[test]
    fn test_fa_cnt() {
        let (reader, _compression) = read_file(Path::new("tests/tgz.fa.gz")).unwrap();
        let fa_records = bio::io::fasta::Reader::new(reader).records();

        for record in fa_records {
            let record = record.unwrap();
            assert_eq!(record.id(), "seqID1");
            assert_eq!(record.seq().to_vec(), b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        }
        
    }
}
*/
pub fn is_fastx(f: String) -> Result<(), String> {
    let ext = vec!["fa", "fas", "fasta", "fastq", "fq", "gz"];
    
    let path = Path::new(&f);
    let f_ext = path.extension().unwrap();
    let f_last = f_ext.to_str().unwrap();

    if ext.contains(&f_last){
        return Ok(());
    } else {
        return Err("Input file is not fasta nor fastq formatted".to_string());
    }
}