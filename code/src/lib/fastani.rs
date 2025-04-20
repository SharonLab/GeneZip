use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;
use polars::datatypes::{DataType, Field, PlSmallStr};
use polars::prelude::{LazyCsvReader, Schema};
use crate::ani_calculator::AniCalculator;

pub struct FastANI {
    results: PathBuf,
}

impl FastANI {
    /*
    Run an ANI calculator using the given parameters
     */
    pub fn run(input: &Path, output: &Path, n_jobs: usize) -> Self {
        match Command::new("fastANI").arg("--rl").arg(input).arg("--ql").arg(input).arg("-t").arg(n_jobs.to_string()).arg("-o").arg(output).spawn() {
            Err(e) => panic!("E: failed to run fastANI due to '{}'", e),
            Ok(mut p) => { p.wait().unwrap(); },
        }

        Self {
            results: output.to_path_buf(),
        }
    }

    /*
    Use a precalculated ANI
     */
    pub fn pre_calculated(results: &Path) -> Self {
        Self {
            results: results.to_path_buf()
        }
    }
}

impl AniCalculator for FastANI {
    fn get_data_frame(&self) -> LazyCsvReader {
        let schema = Schema::from_iter(vec![
            Field::new(PlSmallStr::from("g1"), DataType::String),
            Field::new(PlSmallStr::from("g2"), DataType::String),
            Field::new(PlSmallStr::from("ani"), DataType::Float64),
            Field::new(PlSmallStr::from("bidi_frags"), DataType::Float64),
            Field::new(PlSmallStr::from("total_frags"), DataType::UInt64),
        ]);

        LazyCsvReader::new(self.results.as_path())
            .with_has_header(false)
            .with_separator(b'\t')
            .with_schema(Some(Arc::new(schema)))
    }

    fn results_path(&self) -> &Path { self.results.as_path() }
}