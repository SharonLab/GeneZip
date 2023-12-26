use std::collections::hash_map::IterMut;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::Error;

pub type OutputStream = Box<BufWriter<File>>;

#[derive(Eq, PartialEq, Copy, Clone, Hash)]
pub enum OutputFileType {
    /// Default output file, contig2size2best output
    BaseGz,
    /// Expressive LZ values table
    LzValues,
}

impl ToString for OutputFileType {
    fn to_string(&self) -> String {
        match self {
            OutputFileType::BaseGz => String::from("Base GZ"),
            OutputFileType::LzValues => String::from("LZ Values"),
        }
    }
}

/// A struct to hold output streams along withValues) the associated paths.
pub struct OutputStreams {
    streams: HashMap<OutputFileType, OutputStream>,
    paths: HashMap<OutputFileType, PathBuf>,
}

fn create_file(file_path: &Path) -> Result<OutputStream, Error> {
    File::create(file_path).map(BufWriter::new).map(Box::new)
}

impl OutputStreams {

    pub fn new(stream_paths: &HashMap<OutputFileType, &Path>) -> Result<Self, Error> {
        let mut streams = HashMap::new();
        let mut paths = HashMap::new();

        for (file_type, file_path) in stream_paths {
            streams.insert(*file_type, create_file(file_path)?);
            paths.insert(*file_type, file_path.to_path_buf());
        }

        Ok(Self {
            streams,
            paths
        })
    }

    pub fn stream(&mut self, file_type: &OutputFileType) -> Option<&mut OutputStream> { self.streams.get_mut(file_type) }

    pub fn streams_iter(&mut self) -> IterMut<OutputFileType, OutputStream> {
        self.streams.iter_mut()
    }

    pub fn flush(&mut self) -> std::io::Result<()> {
        for fout in self.streams.values_mut() {
            fout.flush()?;
        }
        Ok(())
    }
}

impl std::fmt::Display for OutputStreams {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        let mut first = true;
        for oft in [OutputFileType::BaseGz, OutputFileType::LzValues] {
            if first {
               first = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{}: {}", oft.to_string(), self.paths.get(&oft).map_or("None".to_string(), |p| p.display().to_string()) )?;
        }

        Ok(())
    }
}