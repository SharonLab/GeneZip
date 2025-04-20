use std::path::Path;
use polars::prelude::LazyCsvReader;

pub trait AniCalculator {
    /*
    Assuming the data frame to contain the columns g1, g2 and ani.
    ANI is assumed to be a value between 0.0 and 100.0.
     */
    fn get_data_frame(&self) -> LazyCsvReader;

    fn results_path(&self) -> &Path;
}