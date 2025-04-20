use std::fmt::{Display, Formatter};
use clap::ValueEnum;

#[derive(Copy, Clone, ValueEnum)]
pub enum AniCalculatorTool {
    FastANI,
    Skani,
}

impl Default for AniCalculatorTool {
    fn default() -> Self {
        Self::FastANI
    }
}

impl Display for AniCalculatorTool {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match self {
            AniCalculatorTool::FastANI => "fastANI",
            AniCalculatorTool::Skani => "skani",
        })
    }
}