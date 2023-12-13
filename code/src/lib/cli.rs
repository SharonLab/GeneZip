//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::path::{Path, PathBuf};
use clap::{Parser, Subcommand};

/*
TODO: edit the help strings
*/

#[derive(Parser)]
#[command(arg_required_else_help(true), author = "Or Leibovich, Yochai Meir, and Itai Sharon", version = option_env!("CARGO_PKG_VERSION").unwrap_or("1.0.0"), about = "GeneZip: LZ78-based classifier of DNA sequences", long_about = None)]
pub struct Cli {

    /// Set an upper limit of threads to use, uses 1 by default. Set to 0 to use as many threads as possible
    #[arg(short = 'j', long = "jobs", value_name = "jobs", default_value_t = 1)]
    jobs: usize,

    /// buffer size used to read files, must be >= 1
    #[arg(short = 'b', long = "buffer", value_name = "buffer", default_value_t = 512)]
    buffer_size: usize,

    /// Verbose run, enable statistics printing
    #[arg(short = 'v', long = "verbose", value_name = "verbose", default_value_t = false)]
    print_statistics: bool,

    #[command(subcommand)]
    commands: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a GeneZip database for repeat prediction
    Build {
        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>[\t<full taxonomy string>]
        #[arg(short = 'i', long = "train", value_name = "training", required = true)]
        training_name2file_file: PathBuf,

        /// Maximum depth allowed for the context tree, must be >= 1. Tested up-to 17
        #[arg(short = 'd', long = "depth", value_name = "depth", default_value_t = 13)]
        max_depth: usize,

        /// Use k-mer frequencies for faster prediction. Value is what k to use, set to 0 to disable. Also requires that taxonomy column is used for model building (see training).
        #[arg(short = 'k', long = "kmer", value_name = "kmer", default_value_t = 4, required = true)]
        kmer_size: usize,

        /// Path to create the GeneZip database in
        #[arg(long = "db", value_name = "db", required = true)]
        db: PathBuf
    },
    /// Classify sequences using GeneZip, use a database created with build
    DBPredict {
        /// A file with the list of fasta files for prediction, format is
        /// <cluster-name>\t<fasta-file>
        #[arg(short = 't', long = "predict", value_name = "predict", required = true)]
        prediction_name2file_file: PathBuf,

        /// Where to print the output file
        #[arg(short = 'o', long = "output", value_name = "output", required = true)]
        out_file: PathBuf,

        /// If given, used as output path file for ANI between sequences and their best-hit.
        /// ANI will be calculated using fastANI.
        /// By default, paths to sequences of representatives will be taken from the database, however, if --train is given, it will override database paths.
        #[arg(long = "ani", value_name = "ani")]
        ani_out_file: Option<PathBuf>,

        /* --------------------------------- */

        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>[\t<full taxonomy string>]
        ///
        /// Note: Will only be used to override paths in the database for ANI prediction.
        #[arg(short = 'i', long = "train", value_name = "training")]
        training_name2file_file: Option<PathBuf>,

        /// Set %G+C distance between model and test genome limit. To disable, set to 100.
        #[arg(long = "gc", value_name = "gc", default_value_t = 10.0)]
        gc_limit: f64,

        /// Path a pre-existing GeneZip database. If given k-mer and max depth will be taken from the database.
        #[arg(long = "db", value_name = "db", required = true)]
        db: PathBuf,

        /// Replace the raw GeneZip value by the following: (GZ(G_1, G_2) + GZ(G_2, G_1)) / ( GZ(G_1, G_1) + GZ(G_2, G_2))
        #[arg(short = 'r', long = "reflect", value_name = "reflect", default_value_t = false, hide = true)]
        reflect: bool,
    },
    /// Classify sequences using GeneZip, build (but don't keep) the database and then predict
    TrainPredict {
        /// A file with the list of fasta files for prediction, format is
        /// <cluster-name>\t<fasta-file>
        #[arg(short = 't', long = "predict", value_name = "predict", required = true)]
        prediction_name2file_file: PathBuf,

        /// Where to print the output file
        #[arg(short = 'o', long = "output", value_name = "output", required = true)]
        out_file: PathBuf,

        /// If given, used as output path file for ANI between sequences and their best-hit.
        /// ANI will be calculated using fastANI.
        #[arg(long = "ani", value_name = "ani")]
        ani_out_file: Option<PathBuf>,

        /* --------------------------------- */

        /// Maximum depth allowed for the context tree, must be >= 1. Tested up-to 17.
        #[arg(short = 'd', long = "depth", value_name = "depth", default_value_t = 13)]
        max_depth: usize,

        /// Set %G+C distance between model and test genome limit. To disable, set to 100.
        #[arg(long = "gc", value_name = "gc", default_value_t = 10.0)]
        gc_limit: f64,

        /// Use k-mer frequencies for better prediction. Value is what k to use, set to 0 to disable. Also requires that kmer cluster column is used for model building.
        #[arg(short = 'k', long = "kmer", value_name = "kmer", default_value_t = 4)]
        kmer_size: usize,

        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>[\t<full taxonomy string>]
        #[arg(short = 'i', long = "train", value_name = "training", required = true)]
        training_name2file_file: PathBuf,

        /// Replace the raw GeneZip value by the following: (GZ(G_1, G_2) + GZ(G_2, G_1)) / ( GZ(G_1, G_1) + GZ(G_2, G_2))
        #[arg(short = 'r', long = "reflect", value_name = "reflect", default_value_t = false, hide = true)]
        reflect: bool,
    },
    /// Print k-mer frequencies to a TSV file
    PrintKmer {
        /// A file with the list of fasta files in the format
        /// <fasta name>\t<fasta-file>
        #[arg(short = 'i', long = "input", value_name = "input", required = true)]
        input: PathBuf,

        /// Output file path
        #[arg(short = 'o', long = "output", value_name = "output", required = true)]
        output: PathBuf,

        /// Value of k for k-mer
        #[arg(short = 'k', value_name = "k", required = true)]
        k: usize,

        /// Ratio, if given, will print value between 0 and 1, else, will print values from 0 to 100.
        #[arg(short = 'r', long = "ratio", value_name = "ratio", default_value_t = false)]
        ratio: bool,
    },
    /// Create k-mer database of reference sequences for later GeneZip prediction call
    #[command(hide = true)]
    BuildKmers {
        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>[\t<full taxonomy string>]
        #[arg(short = 'i', long = "train", value_name = "training", required = true)]
        training_name2file_file: PathBuf,

        /// Use k-mer frequencies for better prediction. Value is what k to use, set to 0 to disable. Also requires that kmer cluster column is used for model building.
        #[arg(short = 'k', long = "kmer", value_name = "kmer", default_value_t = 4, required = true)]
        kmer_size: usize,

        /// Path to create the GeneZip kmers database in
        #[arg(long = "db", value_name = "db", required = true)]
        db: PathBuf
    },
    /// Print k-mer Pearson correlation best-hit table, with similar structure to GeneZip best-hit table.
    #[command(hide = true)]
    KMerPredict {
        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>
        #[arg(short = 'i', long = "train", value_name = "training", required = true)]
        training_name2file_file: PathBuf,

        /// A file with the list of fasta files for prediction, format is
        /// <cluster-name>\t<fasta-file>
        #[arg(short = 't', long = "predict", value_name = "predict", required = true)]
        prediction_name2file_file: PathBuf,

        /// Where to print the output file
        #[arg(short = 'o', long = "output", value_name = "output", required = true)]
        out_file: PathBuf,

        /// Value is what k to use for k-mer frequencies
        #[arg(short = 'k', long = "kmer", value_name = "kmer", default_value_t = 4)]
        kmer_size: usize,
    },
    /// Classify sequences from a metagenomic assembly fasta file
    MetaPredict {
        /// A fasta file
        #[arg(short = 't', long = "predict", value_name = "predict", required = true)]
        prediction_name2file_file: PathBuf,

        /// Where to print the output file
        #[arg(short = 'o', long = "output", value_name = "output", required = true)]
        out_file: PathBuf,

        /* --------------------------------- */

        /// Maximum depth allowed for the context tree, must be >= 1. Tested up-to 17.
        #[arg(short = 'd', long = "depth", value_name = "depth", default_value_t = 13)]
        max_depth: usize,

        /// A file with the list of fasta files for the cluster models in the format
        /// <cluster-name>\t<fasta-file>
        #[arg(short = 'i', long = "train", value_name = "training", required = true)]
        training_name2file_file: PathBuf,

        /// Assume that the file passed gene prediction, sequence ID format must be ID_N where ID is the contig ID and N is a running index. Sequence with identical ID must show up one after the other.
        #[arg(long = "genes", value_name = "genes", default_value_t = false)]
        genes: bool,

        /// Minimal number of genes in a contig for it to be given prediction. Only acive if --genes is set. To disable (by default) set to 0. For example, if set to 3, only contigs that have three or more genes will be evaluated.'
        #[arg(long = "mingenes", value_name = "mingenes", default_value_t = 0)]
        min_genes: usize,
    }
}

struct FeatureSettings {
    training_name2file_file: PathBuf,
    max_depth: usize,
    kmer_size: Option<usize>,
}

impl FeatureSettings {
    fn new(training_name2file_file: &Path, max_depth: usize, kmer_size: Option<usize>) -> Self {
        FeatureSettings {
            training_name2file_file: training_name2file_file.to_path_buf(),
            max_depth,
            kmer_size: match kmer_size {
                Some(k) => if k == 0 { None } else { Some(k) },
                None => None,
            }
        }
    }
}

struct BuildDBSettings {
    db: PathBuf,
    feature_settings: FeatureSettings
}

impl BuildDBSettings {
    fn new(db: &Path, training_name2file_file: &Path, max_depth: usize, kmer_size: usize) -> Self {
        BuildDBSettings {
            db: db.to_path_buf(),
            feature_settings: FeatureSettings::new(training_name2file_file, max_depth, Some(kmer_size)),
        }
    }
}

struct BuildKmerDBSettings {
    db: PathBuf,
    training_name2file_file: PathBuf,
    kmer_size: usize,
}

impl BuildKmerDBSettings {
    fn new(db: &Path, training_name2file_file: &Path, kmer_size: usize) -> Self {
        BuildKmerDBSettings {
            db: db.to_path_buf(),
            training_name2file_file: training_name2file_file.to_path_buf(),
            kmer_size,
        }
    }
}

struct RunSettings {
    version: &'static str,
    jobs: Option<usize>,
    buffer_size: usize,
    print_statistics: bool,
}

impl RunSettings {
    fn new(version: &'static str, jobs: usize, buffer_size: usize, print_statistics: bool) -> Self {
        RunSettings {
            version,
            jobs: if jobs < 2 { None } else { Some(jobs) },
            buffer_size,
            print_statistics
        }
    }
}

struct PredictionSettings {
    prediction_name2file_file: PathBuf,
    out_file: PathBuf,
    ani_out_file: Option<PathBuf>,
    gc_limit: Option<f64>,
    reflect: bool,
}

impl PredictionSettings {
    fn new(prediction_name2file_file: &Path, out_file: &Path, ani_out_file: &Option<PathBuf>, gc_limit: f64, reflect: bool) -> Self {
        PredictionSettings {
            prediction_name2file_file: prediction_name2file_file.to_path_buf(),
            out_file: out_file.to_path_buf(),
            ani_out_file: ani_out_file.clone(),
            gc_limit: if gc_limit == 100.0 {
                None
            } else {
                Some(gc_limit)
            },
            reflect
        }
    }
}


struct PrintKmerSettings {
    input: PathBuf,
    output: PathBuf,
    k: usize,
    ratio: bool,
}

impl PrintKmerSettings {
    fn new(input: &Path, output: &Path, k: usize, ratio: bool) -> Self {
        PrintKmerSettings {
            input: input.to_path_buf(),
            output: output.to_path_buf(),
            k,
            ratio,
        }
    }
}

struct MetaPrediction {
    prediction_name2file_file: PathBuf,
    out_file: PathBuf,
    genes: bool,
    min_genes: usize,
}

impl MetaPrediction {
    fn new(prediction_name2file_file: &Path, out_file: &Path, genes: bool, min_genes: usize) -> Self {
        Self {
            prediction_name2file_file: prediction_name2file_file.to_path_buf(),
            out_file: out_file.to_path_buf(),
            genes,
            min_genes,
        }
    }
}

enum Task {
    BuildDB(BuildDBSettings),
    DBPredict(PathBuf, Option<PathBuf>, PredictionSettings),
    Predict(FeatureSettings, PredictionSettings),
    PrintKmer(PrintKmerSettings),
    BuildKmer(BuildKmerDBSettings),
    KMerPredict(FeatureSettings, PredictionSettings),
    MetaPredict(FeatureSettings, MetaPrediction),
}

impl From<Commands> for Task {
    fn from(commands: Commands) -> Self {
        match commands {
            Commands::Build {training_name2file_file, max_depth, kmer_size, db} => {
                Task::BuildDB(BuildDBSettings::new(&db, &training_name2file_file, max_depth, kmer_size))
            },
            Commands::DBPredict {prediction_name2file_file, out_file, ani_out_file, training_name2file_file, gc_limit, db, reflect} => {
                Task::DBPredict(db,
                                training_name2file_file,
                                PredictionSettings::new(&prediction_name2file_file, &out_file, &ani_out_file, gc_limit, reflect))
            },
            Commands::TrainPredict {prediction_name2file_file, out_file, ani_out_file, max_depth, gc_limit, kmer_size, training_name2file_file, reflect} => {
                Task::Predict(FeatureSettings::new(&training_name2file_file, max_depth, Some(kmer_size)),
                              PredictionSettings::new(&prediction_name2file_file, &out_file, &ani_out_file, gc_limit, reflect))
            },
            Commands::PrintKmer {input, output, k, ratio} => {
                Task::PrintKmer(PrintKmerSettings::new(&input, &output, k, ratio))
            },
            Commands::BuildKmers {training_name2file_file, kmer_size, db} => {
                Task::BuildKmer(BuildKmerDBSettings::new(&db, &training_name2file_file, kmer_size))
            },
            Commands::KMerPredict {training_name2file_file, prediction_name2file_file, out_file, kmer_size} => {
                Task::KMerPredict(FeatureSettings::new(&training_name2file_file, 13, Some(kmer_size)),
                                  PredictionSettings::new(&prediction_name2file_file, &out_file, &None, 100.0, false))
            },
            Commands::MetaPredict {prediction_name2file_file, out_file, max_depth, training_name2file_file, genes, min_genes} => {
                Task::MetaPredict(FeatureSettings::new(&training_name2file_file, max_depth, None),
                                   MetaPrediction::new(&prediction_name2file_file, &out_file, genes, min_genes))
            }
        }
    }
}

pub enum UserTask {
    BuildDB,
    DBPredict,
    Predict,
    PrintKmer,
    BuildKmer,
    KMerPredict,
    MetaPredict,
}

impl From<&Task> for UserTask {
    fn from(value: &Task) -> Self {
        match value {
            Task::Predict(_, _) => { UserTask::Predict },
            Task::DBPredict(_, _, _) => { UserTask::DBPredict },
            Task::BuildDB(_) => { UserTask::BuildDB },
            Task::PrintKmer(_) => { UserTask::PrintKmer }
            Task::BuildKmer(_) => { UserTask::BuildKmer }
            Task::KMerPredict(_, _) => { UserTask::KMerPredict }
            Task::MetaPredict(_, _) => { UserTask::MetaPredict }
        }
    }
}

pub struct Usage {
    run_settings: RunSettings,
    task: Task,
}

impl Usage {
    pub fn new() -> Self {
        let cli = Cli::parse();

        Usage {
            run_settings: RunSettings::new(option_env!("CARGO_PKG_VERSION").unwrap_or("1.0.0"), cli.jobs, cli.buffer_size, cli.print_statistics),
            task: Task::from(cli.commands.unwrap()),
        }
    }

    pub fn get_database_path(&self) -> Option<&Path> {
        match &self.task {
            Task::BuildDB(s) => Some(&s.db),
            Task::DBPredict(s, _, _) => Some(s),
            Task::Predict(_, _) => None,
            Task::PrintKmer(_) => None,
            Task::BuildKmer(s) => Some(&s.db),
            Task::KMerPredict(_, _) => None,
            Task::MetaPredict(_, _) => None,
        }
    }
    pub fn get_task(&self) -> UserTask { UserTask::from(&self.task) }
    pub fn get_training_name2file_file(&self) -> Option<&Path> {
        match &self.task {
            Task::BuildDB(s) => Some(&s.feature_settings.training_name2file_file),
            Task::DBPredict(_, s, _) => s.as_ref().map(|s| s.as_path()),
            Task::Predict(s, _) => Some(&s.training_name2file_file),
            Task::PrintKmer(s) => Some(&s.input),
            Task::BuildKmer(s) => Some(&s.training_name2file_file),
            Task::KMerPredict(s, _) => Some(&s.training_name2file_file),
            Task::MetaPredict(s, _) => Some(&s.training_name2file_file),
        }
    }
    pub fn get_prediction_name2file_file(&self) -> Option<&Path> {
        match &self.task {
            Task::BuildDB(_) => None,
            Task::DBPredict(_, _, s) => Some(&s.prediction_name2file_file),
            Task::Predict(_, s) => Some(&s.prediction_name2file_file),
            Task::PrintKmer(_) => None,
            Task::BuildKmer(_) => None,
            Task::KMerPredict(_, s) => Some(&s.prediction_name2file_file),
            Task::MetaPredict(_, s) => Some(&s.prediction_name2file_file),
        }
    }
    pub fn get_out_file(&self) -> Option<&Path> {
        match &self.task {
            Task::BuildDB(_) => None,
            Task::DBPredict(_, _, s) => Some(&s.out_file),
            Task::Predict(_, s) => Some(&s.out_file),
            Task::PrintKmer(s) => Some(&s.output),
            Task::BuildKmer(_) => None,
            Task::KMerPredict(_, s) => Some(&s.out_file),
            Task::MetaPredict(_, s) => Some(&s.out_file),
        }
    }
    pub fn get_max_depth(&self) -> Option<usize> {
            match &self.task {
                Task::BuildDB(s) => Some(s.feature_settings.max_depth),
                Task::DBPredict(_, _, _) => None,
                Task::Predict(s, _) => Some(s.max_depth),
                Task::PrintKmer(_) => None,
                Task::BuildKmer(_) => None,
                Task::KMerPredict(_, _) => None,
                Task::MetaPredict(s, _) => Some(s.max_depth),
            }
    }
    pub fn get_version(&self) -> &str { self.run_settings.version }
    pub fn get_jobs(&self) -> Option<usize> { self.run_settings.jobs }
    pub fn get_buffer_size(&self) -> usize { self.run_settings.buffer_size }
    pub fn get_print_statistics(&self) -> bool { self.run_settings.print_statistics }
    pub fn get_gc_limit(&self) -> Option<f64> {
        match &self.task {
            Task::BuildDB(_) => None,
            Task::DBPredict(_, _, s) => s.gc_limit,
            Task::Predict(_, s) => s.gc_limit,
            Task::PrintKmer(_) => None,
            Task::BuildKmer(_) => None,
            Task::KMerPredict(_, _) => None,
            Task::MetaPredict(_, _) => None,
        }
    }
    pub fn get_ani_out_file(&self) -> Option<&Path> {
        match &self.task {
            Task::BuildDB(_) => None,
            Task::DBPredict(_, _, s) => match &s.ani_out_file {
                Some(p) => Some(p),
                None => None,
            },
            Task::Predict(_, s) => match &s.ani_out_file {
                Some(p) => Some(p),
                None => None,
            },
            Task::PrintKmer(_) => None,
            Task::BuildKmer(_) => None,
            Task::KMerPredict(_, _) => None,
            Task::MetaPredict(_, _) => None,
        }
    }
    pub fn get_kmer_size(&self) -> Option<usize> {
        match &self.task {
            Task::BuildDB(s) => s.feature_settings.kmer_size,
            Task::Predict(s, _) => s.kmer_size,
            Task::DBPredict(_, _, _) => None,
            Task::PrintKmer(s) => Some(s.k),
            Task::BuildKmer(s) => Some(s.kmer_size),
            Task::KMerPredict(s, _) => s.kmer_size,
            Task::MetaPredict(_, _) => None,
        }
    }

    pub fn get_ratio(&self) -> bool {
        match &self.task {
            Task::PrintKmer(s) => s.ratio,
            _ => false,
        }
    }

    pub fn get_reflect(&self) -> bool {
        match &self.task {
            Task::Predict(_, s) => s.reflect,
            Task::DBPredict(_, _, s) => s.reflect,
            _ => false,
        }
    }

    pub fn get_genes(&self) -> bool {
        match &self.task {
            Task::MetaPredict(_, s) => s.genes,
            _ => false,
        }
    }

    pub fn get_min_genes(&self) -> usize {
        match &self.task {
            Task::MetaPredict(_, s) => s.min_genes,
            _ => 0,
        }
    }
}
