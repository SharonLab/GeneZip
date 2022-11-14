use std::path::{Path, PathBuf};
use clap::{Arg, App, value_t};


pub struct Usage {
    max_depth: usize,
    training_name2file_file: PathBuf,
    prediction_name2file_file: PathBuf,
    out_file: PathBuf,
    version: &'static str,
    jobs: Option<usize>,
    buffer_size: usize,
}

impl Usage {
    pub fn new() -> Self {
        let max_depth_was_tested_upto = 17_usize;
        let default_buffer_size = "512";
        let default_max_depth = "13";
        let default_num_threads = "1";
        let version = option_env!("CARGO_PKG_VERSION").unwrap_or("1.0.0"); // TODO: use cargo version

        let depth_help = format!("maximum depth allowed for the context tree, must be >= 1 (default: {}). Tested up-to {}.", default_max_depth, max_depth_was_tested_upto);
        let buffer_size_help = format!("buffer size used to read files, must be >= 1 (default: {})", default_buffer_size);
        let jobs_help = format!("Set an upper limit of threads to use, uses 1 by default. Set to 0 to use as many threads as possible.");

        let pre_matches = App::new("GeneZip")
            .version(version)
            .name("GeneZip")
            .arg(Arg::with_name("training")
                .takes_value(true)
                .short('i')
                .long("training")
                .help("a file with the list of fasta files for the cluster models in the format\n<cluster-name>\t<fasta-file>")
                .required(true))
            .arg(Arg::with_name("predict")
                .name("predict")
                .short('t')
                .takes_value(true)
                .required(true)
                .help("a file with the list of fasta files for prediction, format is\n<cluster-name>\t<fasta-file>"))
            .arg(Arg::with_name("output")
                .name("output")
                .short('o')
                .takes_value(true)
                .required(true)
                .help("name of the output file"))
            .arg(Arg::with_name("depth")
                .name("depth")
                .short('d')
                .takes_value(true)
                .default_value(default_max_depth)
                .help(depth_help.as_str()))
            .arg(Arg::with_name("buffer")
                .name("buffer")
                .takes_value(true)
                .default_value(default_buffer_size)
                .help(buffer_size_help.as_str()))
            .arg(Arg::with_name("jobs")
                .long("jobs")
                .short('j')
                .default_value(default_num_threads)
                .help(jobs_help.as_str()));

        let matches = pre_matches.get_matches();

        let buffer: usize = value_t!(matches.value_of("buffer"), usize).unwrap();
        if buffer < 1 {
            panic!("E: buffer size (--buffer) should not be below 1, but given value is {}", buffer);
        }
        let depth: usize = value_t!(matches.value_of("depth"), usize).unwrap();
        if depth < 1 {
            panic!("E: depth (-d) should not be below 1, but given value is {}", depth);
        }
        let training: PathBuf = Path::new(matches.value_of("training").unwrap()).to_path_buf();
        let predict: PathBuf = Path::new(matches.value_of("predict").unwrap()).to_path_buf();
        let output: PathBuf = Path::new(matches.value_of("output").unwrap()).to_path_buf();

        let raw_jobs: usize = value_t!(matches.value_of("jobs"), usize).unwrap();

        if depth > max_depth_was_tested_upto {
            eprintln!("W: max depth is set to {}, however GeneZip has been tested up-to with max depth up-to {}.", depth, max_depth_was_tested_upto)
        }

        Usage {
            max_depth: depth,
            training_name2file_file: training,
            prediction_name2file_file: predict,
            out_file: output,
            version,
            jobs: match raw_jobs {
                0 => None,
                _ => Some(raw_jobs)
            },
            buffer_size: buffer,
        }
    }

    pub fn get_training_name2file_file(&self) -> &Path { &self.training_name2file_file }
    pub fn get_prediction_name2file_file(&self) -> &Path { &self.prediction_name2file_file }
    pub fn get_out_file(&self) -> &Path { &self.out_file }
    pub fn get_max_depth(&self) -> usize { self.max_depth }
    pub fn get_version(&self) -> &str { self.version }
    pub fn get_jobs(&self) -> Option<usize> { self.jobs }
    pub fn get_buffer_size(&self) -> usize { self.buffer_size }
}
