/*
Prepare reference database for GeneZip

1. Collect N random samples from each species
2. For each genus, calculate ANI between all references
    2.1. if any of them are too similar (distance <= 5), merge them
*/
use GeneZipLib::skani::Skani;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::collections::{HashSet, HashMap};
use polars::prelude::*;
use polars::prelude::Expr::Column;
use rand::rngs::StdRng;
use rand::prelude::IndexedRandom;

use clap::Parser;
use rand::SeedableRng;
use GeneZipLib::graph::Graph;
use GeneZipLib::static_graph::StaticGraph;
use GeneZipLib::fastani::FastANI;
use GeneZipLib::ani_calculator::AniCalculator;

use std::option::Option;
use GeneZipLib::ani_calculator_tool::AniCalculatorTool;
use GeneZipLib::samples_file_reader::{Sample, SampleSource};
use GeneZipLib::taxonomy::{TaxonomicRank, Taxonomy};

type TaxonomyName = std::string::String;
// type MagName = std::string::String;

// Collect the full sample provided by the user
// Map taxonomies to the order in which they were found, this can be used to keep consistent order
fn collect_samples(input_database: &Path) -> (HashMap<Taxonomy, Vec<Sample>>, HashMap<Taxonomy, usize>) {
    let mut mapping = HashMap::new();
    let mut order = HashMap::new();

    for sample in SampleSource::new(input_database, true).into_iter() {
        match sample {
            Err(e) => panic!("E: failed to collect samples from '{}' due to: {}", input_database.display(), e),
            Ok(sample) => {
                let taxonomy = sample.get_taxonomy().clone();
                match taxonomy {
                    Some(mut taxonomy) => {
                        if !taxonomy.has_rank(&TaxonomicRank::Species) {
                            taxonomy.fill_rank(&TaxonomicRank::Species);
                        }
                        mapping.entry(taxonomy.clone()).or_insert(Vec::new()).push(sample);
                        order.entry(taxonomy).or_insert(mapping.len() - 1);
                    },
                    None => panic!("E: invalid sample has no taxonomy"),
                }
            },
        }
    }

    (mapping, order)
}

fn taxa2order_into_ordered_taxa(taxa2order: &HashMap<Taxonomy, usize>) -> Vec<&Taxonomy> {
    let order2taxa = taxa2order.iter().map(|(k, v)| (*v, k)).collect::<HashMap<usize, &Taxonomy>>();
    let mut order = order2taxa.keys().cloned().collect::<Vec<usize>>();
    order.sort();
    order.iter().map(|o| order2taxa[o]).collect()
}

// Randomly pick N (number_of_representatives_per_species) representatives from each taxonomy.
fn pick_representatives(mapping: HashMap<Taxonomy, Vec<Sample>>, taxa2order: &HashMap<Taxonomy, usize>, number_of_representatives_per_species: usize, random_state: &mut StdRng) -> HashMap<Taxonomy, Vec<Sample>> {
    // Make sure we sample in a consistent order to achieve repeatable runs
    let ordered_taxa = taxa2order_into_ordered_taxa(taxa2order);
    
    // Now collect
    ordered_taxa.iter()
        .map(|&tax| (tax.clone(), mapping[tax].choose_multiple(random_state, number_of_representatives_per_species).cloned().collect::<Vec<Sample>>()))
        .collect()
}

fn apply_func_on_str_series(str_series: polars::frame::column::Column, func: &dyn Fn (&str) -> String) -> PolarsResult<Option<polars::frame::column::Column>> {
    let results = str_series.
        str()?
        .into_iter()
        .map(|s: Option<&str> | s.map(func))
        .collect::<StringChunked>();
    Ok(Some(results.into_column()))
}

// fn merge_by_ani(ani_results_path: &dyn AniCalculator, path_str2taxonomy: HashMap<String, Taxonomy>, ani_cutoff: f64) -> PolarsResult<Vec<HashSet<Taxonomy>>> {
fn merge_by_ani(ani_results_path: &dyn AniCalculator, mag_path2taxonomy: &HashMap<&Path, &Taxonomy>, ani_cutoff: f64) -> PolarsResult<Vec<HashSet<Taxonomy>>> {
    let g1 = Column(PlSmallStr::from("g1"));
    let g2 = Column(PlSmallStr::from("g2"));
    let ani = Column(PlSmallStr::from("ani"));

    let mag_path2taxonomy = mag_path2taxonomy.iter().map(|(&k, &v)| (k.to_path_buf(), v.clone())).collect::<HashMap<PathBuf, Taxonomy>>();
    
    let ani_df = {
        let fs = move |s: &str| match mag_path2taxonomy.get(PathBuf::from(s).as_path()) {
            Some(taxonomy) => taxonomy.to_string(),
            None => panic!("E: '{}' doesn't represent a path in mag_path2taxonomy, can't recover", s),
        };
        let f = move |c| apply_func_on_str_series(c, &fs);
        ani_results_path.get_data_frame().finish()?
            // .filter(g1.clone().lt_eq(g2.clone())) // This leads to lost taxa (make sense)
            .filter((lit(100.0) - ani).lt_eq(ani_cutoff))
            .with_column(g1.clone().map(f.clone(), GetOutput::same_type()).alias(PlSmallStr::from("g1")))
            .with_column(g2.clone().map(f, GetOutput::same_type()).alias(PlSmallStr::from("g2")))
            .select([g1, g2])
            .collect()
    }?;

    let mut graph: StaticGraph<TaxonomyName, u8, u8> = StaticGraph::new(None, 0);
    let mut row = ani_df.get_row(0)?;
    for i in 0..ani_df.height() {
        ani_df.get_row_amortized(i, &mut row)?;
        let g1_data = row.0.get(0).unwrap_or_else(|| panic!("E: Row {} is missing g1 data", i)).to_string();
        let g2_data = row.0.get(1).unwrap_or_else(|| panic!("E: Row {} is missing g2 data", i)).to_string();
        graph.add_edge(&g1_data, &g2_data, &None, |_, _| None);
    }

    Ok(graph.connected_components().iter().map(|s| s.iter().map(|&i| Taxonomy::from(i.as_str())).collect()).collect())
}

fn get_genus_cluster_list_path(work_folder: &Path, genus_number: usize) -> PathBuf {
    PathBuf::from(work_folder).join(format!("{:?}.genus_cluster.list", genus_number))
}

fn get_genus_cluster_ani_path(work_folder: &Path, genus_number: usize) -> PathBuf {
    PathBuf::from(work_folder).join(format!("{:?}.genus_cluster.ani", genus_number))
}

// Transform (full) Taxonomy2Samples mapping into Sample2Genus mapping, discarding the rank of species.
// Also returning a set of found genera.
// This function create a file for each genus that contain all its samples.
fn get_mag_path2taxonomy_and_found_genera<'a>(mapping: &'a HashMap<Taxonomy, Vec<Sample>>,
                                              taxa2order: &'a HashMap<Taxonomy, usize>,
                                              work_folder: &Path) -> std::io::Result<HashMap<&'a Path, &'a Taxonomy>> {
    let mut mag_path2taxonomy = HashMap::new();

    // Create list files for ani runs
    let mut found_genera = HashSet::new();

    let mut genus2handler: HashMap<Taxonomy, BufWriter<File>> = HashMap::new();
    let genus2order = &get_genus_order(taxa2order);
    
    // Using ordered taxa will lead to consistent genus list files
    let ordered_taxa = taxa2order_into_ordered_taxa(taxa2order);
    
    for taxonomy in ordered_taxa {
        let mags_vec = &mapping[&taxonomy];
        let genus = taxonomy.limit2rank(&TaxonomicRank::Genus);
        
        found_genera.insert(genus.clone());
        let genus_output_fp = get_genus_cluster_list_path(work_folder, *genus2order.get(&genus).unwrap_or_else(|| panic!("E: genus '{}' have no order!", genus)));
        if !genus2handler.contains_key(&genus) {
            genus2handler.insert(genus.clone(),
                                 BufWriter::new(File::create_new(genus_output_fp.clone())?));
        }
        assert!(genus_output_fp.exists());
        
        for sample in mags_vec {
            mag_path2taxonomy.insert(sample.get_path(), taxonomy);
            writeln!(genus2handler.get_mut(&genus).expect("E: unreachable"), "{}", sample.get_path().display())?;
        }
    }
    assert_eq!(found_genera.len(), genus2order.len());

    for fh in genus2handler.values_mut() {
        fh.flush()?;
    }

    Ok(mag_path2taxonomy)
}

fn get_ani_calculators(taxa2order: &HashMap<Taxonomy, usize>, work_folder: &Path, n_jobs: usize, ani_calculator: AniCalculatorTool) -> Vec<Box<dyn AniCalculator>> {
    taxa2order.values()
        .filter_map(|i| {
            let list_path = get_genus_cluster_list_path(work_folder, *i);
            let ani_results_path = get_genus_cluster_ani_path(work_folder, *i);
            let ani: Option<Box<dyn AniCalculator>> = if ! ani_results_path.exists() {
                if list_path.exists() {
                    Some(match ani_calculator {
                        AniCalculatorTool::FastANI => Box::new(FastANI::run(get_genus_cluster_list_path(work_folder, *i).as_path(),
                                                                            ani_results_path.as_path(),
                                                                            n_jobs)),
                        AniCalculatorTool::Skani => Box::new(Skani::run(get_genus_cluster_list_path(work_folder, *i).as_path(),
                                                                        ani_results_path.as_path(),
                                                                        n_jobs)),
                    })
                } else {
                    None
                }
            } else {
                Some(match ani_calculator {
                    AniCalculatorTool::FastANI => Box::new(FastANI::pre_calculated(ani_results_path.as_path())),
                    AniCalculatorTool::Skani => Box::new(Skani::pre_calculated(ani_results_path.as_path())),
                })
            };
            ani
        })
        .collect()
}

// Merge taxonomies by ANI values.
fn merge_by_genus(mapping: HashMap<Taxonomy, Vec<Sample>>,
                  taxa2order: &HashMap<Taxonomy, usize>,
                  merge_distance: f64,
                  work_folder: &Path,
                  n_jobs: usize,
                  ani_calculator: AniCalculatorTool) -> std::io::Result<HashMap<Vec<Taxonomy>, HashSet<Sample>>> {
    let mag_path2taxonomy= get_mag_path2taxonomy_and_found_genera(&mapping, taxa2order, work_folder)?;
    

    // Run ANI calculator
    let ani_calculators = get_ani_calculators(taxa2order, work_folder, n_jobs, ani_calculator);

    // Extract merges
    let merges: Vec<HashSet<Taxonomy>> = ani_calculators.iter()
        .flat_map(|ani| merge_by_ani(ani.as_ref(), &mag_path2taxonomy, merge_distance)
            .unwrap_or_else(|e| panic!("E: failed to extract data from '{}' because '{:?}'", ani.results_path().display(), e)))
        .collect();

    // Create results
    let mut results: HashMap<Vec<Taxonomy>, HashSet<Sample>> = HashMap::new();
    let mut merged_taxa: HashSet<Taxonomy> = HashSet::new();
    for mer in merges.iter().map(|m| Vec::from_iter(m.iter().cloned())) {
        results.insert(mer.clone(), HashSet::new()); // mer is turned into a tuple on the python version
        for m in mer.iter() {
            results.get_mut(&mer).expect("E: we did insert the empty set, so it should wrong.").extend(mapping[m].iter().cloned());
        }
        merged_taxa.extend(mer.iter().cloned());
    }

    for (m, k) in mapping {
        if ! merged_taxa.contains(&m) {
            results.insert(vec![m.clone()], HashSet::from_iter(k)); // m was in a tuple
        }
    }

    Ok(results)
}

#[allow(dead_code)]
fn drop_no_species(mapping: &mut HashMap<Taxonomy, Vec<Sample>>) {
    let mut to_remove = Vec::new();
    for taxonomy in mapping.keys() {
        let empty = match taxonomy.get_taxa(&TaxonomicRank::Species) {
            None => true,
            Some(s) => s.is_empty(),
        };
        if empty {
            to_remove.push(taxonomy.clone());
        }
    }

    for tax in to_remove {
        mapping.remove(&tax);
    }
}

fn create_taxa2cluster(merged_taxa2mags: &HashMap<Vec<Taxonomy>, HashSet<Sample>>, work_folder: &Path) -> std::io::Result<()> {
    let path = PathBuf::from(work_folder).join("taxa2cluster.tsv");
    let mut stream = BufWriter::new(File::create_new(path.as_path())?);
    for (i, cluster) in merged_taxa2mags.keys().enumerate() {
        for member in cluster {
            writeln!(stream, "{}\t{}", member, i)?;
        }
    }
    stream.flush()
}

fn create_representative2cluster(merged_taxa2mags: &HashMap<Vec<Taxonomy>, HashSet<Sample>>, work_folder: &Path) -> std::io::Result<()> {
    let path = PathBuf::from(work_folder).join("representative2cluster.tsv");
    let mut stream = BufWriter::new(File::create_new(path.as_path())?);
    for (i, cluster) in merged_taxa2mags.values().enumerate() {
        for member in cluster {
            writeln!(stream, "{}\t{}", member.get_path().display(), i)?;
        }
    }
    stream.flush()
}

fn create_training_file(merged_taxa2mags: &HashMap<Vec<Taxonomy>, HashSet<Sample>>, work_folder: &Path) -> std::io::Result<()> {
    let path = PathBuf::from(work_folder).join("training.tsv");
    let mut stream = BufWriter::new(File::create_new(path.as_path())?);
    for cluster in merged_taxa2mags.values() {
        for member in cluster {
            writeln!(stream, "{}\t{}", member.get_name(), member.get_path().display())?;
        }
    }
    stream.flush()
}

fn get_genus_order(taxa2order: &HashMap<Taxonomy, usize>) -> HashMap<Taxonomy, usize> {
    // Using ordered taxa will lead to consistent genus order, and we will place lower first to pack the numbers in the names of files
    let ordered_taxa = taxa2order_into_ordered_taxa(taxa2order);
    let mut genus_order = HashMap::new();
    for taxonomy in ordered_taxa {
        let order = taxa2order.get(taxonomy).unwrap();
        let entry = genus_order.entry(taxonomy.limit2rank(&TaxonomicRank::Genus)).or_insert(*order);
        // *entry = (*entry).min(*order); // it should already be the case
        assert!(*entry <= *order);
    }

    genus_order
}

fn run(input_database: &Path, merge_distance: f64, work_folder: &Path, n_jobs: usize, sample_size: usize, random_state: &mut StdRng, ani_calculator: AniCalculatorTool) -> std::io::Result<()> {
    if !work_folder.exists() {
        create_dir_all(work_folder)?;
    }

    let (taxa2mags, taxa2order) = collect_samples(input_database);
    assert_eq!(taxa2mags.len(), taxa2order.len());
    // drop_no_species(&mut taxa2mags);
    let taxa2mags = pick_representatives(taxa2mags, &taxa2order, sample_size, random_state);

    // Debug missing list files:
    // for (taxa, order) in taxa2order.iter() {
    //     eprintln!("T2O:\t{}\t{}", order, taxa);
    // }
    // for (taxa, samples) in taxa2mags.iter() {
    //     eprintln!("T2M:\t{}\t{:?}", taxa, samples);
    // }
    
    assert_eq!(taxa2mags.len(), taxa2order.len());
    let merged_taxa2mags = merge_by_genus(taxa2mags, &taxa2order, merge_distance, work_folder, n_jobs, ani_calculator)?;

    // Create output files
    create_taxa2cluster(&merged_taxa2mags, work_folder)?;
    create_representative2cluster(&merged_taxa2mags, work_folder)?;
    create_training_file(&merged_taxa2mags, work_folder)
}

#[derive(Parser)]
#[command(arg_required_else_help(true), author = "Or Leibovich, Yochai Meir, and Itai Sharon", version = option_env!("CARGO_PKG_VERSION").unwrap_or("1.0.0"), about = "Training database builder for GeneZip", long_about = None)]
pub struct DBBuilderCLI {

    /// Set an upper limit of threads to use, uses 1 by default. Set to 0 to use as many threads as possible
    #[arg(short = 'j', long = "jobs", value_name = "jobs", default_value_t = 1)]
    jobs: usize,

    /// Sample size per species
    #[arg(short = 's', long = "sample", value_name = "sample", required = true)]
    sample: usize,

    /// Output folder
    #[arg(short = 'o', long = "output", value_name = "output", required = true)]
    output: PathBuf,

    /// Merge distance
    #[arg(short = 'm', long = "merge", value_name = "merger", required = true)]
    merge: f64,

    /// Input database
    #[arg(short = 'i', long = "input", value_name = "input", required = true)]
    input: PathBuf,

    /// Random state seed
    #[arg(short = 'r', long = "rss", value_name = "rss", default_value_t = 1)]
    rss: u64,
    
    /// Choose ANI calculation tool
    #[arg(short = 'a', long = "ani", value_name = "ani")]
    ani: AniCalculatorTool,
}
fn interface() {
    let cli = DBBuilderCLI::parse();
    match run(cli.input.as_path(), cli.merge, cli.output.as_path(), cli.jobs, cli.sample, &mut StdRng::seed_from_u64(cli.rss), cli.ani) {
        Ok(_) => eprintln!("Done"),
        Err(e) => eprintln!("E: failed to run to conclusion, got the following error: {}", e),
    }
}


fn main() {
    interface()
}
