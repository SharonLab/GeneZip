use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;
use polars::datatypes::{DataType, Field, PlSmallStr};
use polars::prelude::{LazyCsvReader, Schema};
use crate::ani_calculator::AniCalculator;

pub struct Skani {
    results: PathBuf,
}

impl Skani {
    /*
    Run an ANI calculator using the given parameters
     */
    pub fn run(input: &Path, output: &Path, n_jobs: usize) -> Self {
        match Command::new("skani").arg("dist").arg("-t").arg(n_jobs.to_string()).arg("-o").arg(output).arg("--ql").arg(input).arg("--rl").arg(input).spawn() {
            Err(e) => panic!("E: failed to run skani due to '{}'", e),
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

impl AniCalculator for Skani {
    fn get_data_frame(&self) -> LazyCsvReader {
        let schema = Schema::from_iter(vec![
            // Real names
            // Field::new(PlSmallStr::from("Ref_file"), DataType::String),
            // Field::new(PlSmallStr::from("Query_file"), DataType::String),
            // Field::new(PlSmallStr::from("ANI"), DataType::Float64),
            // Field::new(PlSmallStr::from("Align_fraction_ref"), DataType::Float64),
            // Field::new(PlSmallStr::from("Align_fraction_query"), DataType::Float64),
            // Field::new(PlSmallStr::from("Ref_name"), DataType::String),
            // Field::new(PlSmallStr::from("Query_name"), DataType::String),
            
            // Names consistent with ANI calculator documentation module.
            Field::new(PlSmallStr::from("g1"), DataType::String),
            Field::new(PlSmallStr::from("g2"), DataType::String),
            Field::new(PlSmallStr::from("ani"), DataType::Float64),
            Field::new(PlSmallStr::from("Align_fraction_ref"), DataType::Float64),
            Field::new(PlSmallStr::from("Align_fraction_query"), DataType::Float64),
            Field::new(PlSmallStr::from("Ref_name"), DataType::String),
            Field::new(PlSmallStr::from("Query_name"), DataType::String),
        ]);

        LazyCsvReader::new(self.results.as_path())
            .with_has_header(true)
            // .with_skip_lines(1) // This will skip the first line after the header, not the header.
            .with_separator(b'\t')
            .with_schema(Some(Arc::new(schema)))
    }

    fn results_path(&self) -> &Path { self.results.as_path() }
}