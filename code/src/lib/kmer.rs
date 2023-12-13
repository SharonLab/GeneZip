use std::collections::HashMap;
use ndarray::prelude::*;
use lazy_static::*;
use rayon::prelude::*;
use std::path::Path;
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;

lazy_static! {
    pub static ref ALPHABET_MAP: HashMap<u8, u8> = {
        let mut map = HashMap::new();
        map.insert(b'A', b'T');
        map.insert(b'T', b'A');
        map.insert(b'C', b'G');
        map.insert(b'G', b'C');
        map
    };
}

lazy_static! {
    static ref ALPHABET_VALUES: HashMap<u8, usize> = {
        let mut map = HashMap::new();
        map.insert(b'A', 0b00);
        map.insert(b'T', 0b11);
        map.insert(b'G', 0b10);
        map.insert(b'C', 0b01);
        map
    };
}

lazy_static! {
    pub static ref REV_ALPHABET_VALUES: HashMap<usize, u8> = {
        let mut map = HashMap::new();
        map.insert(0b00, b'A');
        map.insert(0b11, b'T');
        map.insert(0b10, b'G');
        map.insert(0b01, b'C');
        map
    };
}

lazy_static! {
    static ref KMASK: HashMap<usize, usize> = {
        let mut map = HashMap::new();
        for i in 1..=32 {
            map.insert(i, create_bits_mask(i));
        }
        map
    };
}

/// Creates a bit mask to one pick and index K characters.
///
/// Built for 4 characters.
///
/// # Example
///
/// let results = create_bits_mask(2);
/// assert_eq!(0b1111, results);
fn create_bits_mask(k: usize) -> usize {
    let mut mask = 0;
    for _ in 0..k {
        mask = (mask << 2) + 0b11;
    }
    mask
}

/// Calculates the index of reverse complement word.
///
/// Assuming all k-mers are sitting in a vector, each has one index.
/// Meaning, each word, when translated will turn into the index of the reverse-complement word.
///
/// Takes the index of current word and k.
/// Assuming letters can be translated with XOR.
///
/// # Example
///
/// let reverse = complement_index(0b1000, 2);
/// assert_eq!(reverse, 0b1101);
/// let reverse = complement_index(0b0110, 2);
/// assert_eq!(reverse, 0b0110);
fn complement_index(word: usize, k: usize) -> usize {
    let mask: usize = KMASK[&k];
    let nucleotide_mask: usize = KMASK[&1]; // one nucleotide
    let xor_word = word ^ mask;
    let mut new_index: usize = 0;
    for j in (1..=k).rev() {
        new_index = (new_index << 2) | (nucleotide_mask & (xor_word >> ((k - j) * 2)));
    }
    new_index
}

/// Returns a vector of the text-side of the genomic profiles, meaning, the kept kmers after
/// the vector minimization.
pub fn create_printed_vector(k: usize, alphabet_values: &HashMap<u8, u8>, rev_alphabet_values: &HashMap<usize, u8>) -> Vec<String> {
    let mut results = Vec::new();
    let nucleotide_mask = create_bits_mask(1);
    for i in 0..alphabet_values.len().pow(k as u32) {

        // Now find it counter part and remove from the vector
        if complement_index(i, k) >= i { // We need to keep this string
            // Translate the index i into the sequence that would have been counted here.
            let mut sequence = String::new();
            let mut masked_i = i;
            for _ in 0..k {
                let key = masked_i & nucleotide_mask;
                sequence.insert(0, rev_alphabet_values[&key] as char);
                masked_i >>= 2;
            }
            results.push(sequence.clone());
        }
    }

    results
}


/// Helper function to get_fasta_vector, allow abstraction of BufRead to open text file, gz file etc...
fn get_fasta_vector<I>(k: usize, fasta: I, alphabet_map: &HashMap<u8, u8>, alphabet_values: &HashMap<u8, usize>) -> Result<Array1<f64>, String>
    where I: IntoIterator<Item=u8> {
    let mut vector = Array::zeros(alphabet_values.len().pow(k as u32));
    let mask = KMASK[&k];
    let mut index: usize = 0; // left to right word index
    let mut findex: usize = 0; // flipped word index
    let mut found = 0; // how many nucleotides read so far, we own't start counting before k

    for nucleotide in fasta {
        if alphabet_map.keys().any(|x| *x == nucleotide) {
            // Handle first vector
            index = (index << 2) | (alphabet_values[&nucleotide]);
            index &= mask;

            // Handle flipped-words vector
            findex = (findex >> 2) | ((alphabet_values[&alphabet_map[&nucleotide]]) << (2 * (k - 1)));


            if found < k {
                found += 1;
            }
            if k <= found {
                if index < findex {
                    vector[index] += 1.0; // Count the found word
                } else {
                    vector[findex] += 1.0; // Count the flipped word
                }
            }
        } else if nucleotide == b'N' {
            found = 0;
        }
    }

    Ok(vector)
}

/// Minimizes the given vector by dropping the reverse-complement words.
/// For example k=4, AATC <=> GATT, so only one of them will be kept.
///
/// # Arguments
///
/// * `vector` - A vector to minimize
/// * `k` - How long should each kmer should be be.
///
fn minimize_vector(vector: &mut Array1<f64>, k: usize) -> Array1<f64> {
    let mut new_vector_size = 0;
    for i in 0..vector.shape()[0] {
        if vector[i] >= 0.0 {
            // Now find it counter part and remove from the vector
            let new_index = complement_index(i, k);
            if i != new_index {
                vector[new_index] = -1.0;
            }
            new_vector_size += 1;
        }
    }

    // Count how many items equal-to or above 0, using this count as the size of the new array.
    // let new_vector_size = vector.map(|x| if *x >= 0.0 {1} else {0}).scalar_sum();
    let mut results_array = Array::zeros(new_vector_size);
    let mut i: usize = 0;
    for elem in vector {
        if *elem >= 0.0 {
            results_array[i] = *elem;
            i += 1;
        }
    }

    results_array
}

/// Normalizes the vector values
fn normalize_vector(vector: &mut Array1<f64>) -> &mut Array1<f64> {
    let total_kmer_values = vector.sum();
    *vector /= total_kmer_values;
    vector
}

/// Creates a normalized genomic profile using k kmer.
/// Returns a tuple of fasta path and the vector to allow multithreaded mapping into an HashMap.
///
/// # Arguments
///
/// * `k` - length of kmer.
/// * `fasta_path` - Path to a fasta file.
///
pub fn create_normalized_profile<I, N: ?Sized>(k: usize, fasta: I, id: &N) -> (&N, Result<Array1<f64>, String>)
    where I: IntoIterator<Item=u8> + Clone + Sync {
    if k > usize::MAX / 2 {  // We use 2 bits per nucleotide
        panic!("K-mer value is too high, value must be <= {} on your system", usize::MAX / 2);
    }

    let mut vector = match &mut get_fasta_vector(k, fasta, &ALPHABET_MAP, &ALPHABET_VALUES) {
        Ok(v) => minimize_vector(v, k),
        Err(e) => return (id, Err(e.to_string())),
    };
    normalize_vector(&mut vector);
    (id, Ok(vector))
}

/// Translates the given genomes into k genomic-profile HashMap using up-to ncores but not less then 1.
/// Errors will be returned per-given fasta to avoid stopping the whole computation over one bad file.
///
/// Uses as a choke-point to control the execution of internal functions.
///
/// # Arguments
/// * `genomes` - which genomes to take
/// * `k` - the k in kmer
#[allow(dead_code)]
pub fn generate_kmers<'a>(genomes: &Vec<&'a Path>, k: usize, use_multithreading: bool, buffer_size: usize) -> HashMap<&'a Path, Result<Array1<f64>, String>> {
    let mut result = HashMap::new();
    if use_multithreading {
        result.par_extend(genomes.into_par_iter().map(|&genome| create_normalized_profile(k, FastaNucltudiesIterator::new(genome, buffer_size), genome)));
    } else {
        for genome in genomes {
            let (_, normalized_profile) = create_normalized_profile(k, FastaNucltudiesIterator::new(genome, buffer_size), genome);
            result.insert(*genome, normalized_profile);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;
    use crate::kmer::{ALPHABET_MAP, ALPHABET_VALUES, complement_index, create_normalized_profile, create_printed_vector, get_fasta_vector, minimize_vector, REV_ALPHABET_VALUES};

    #[test]
    fn one_mer() {
        let dna_test_path = PathBuf::from("../tests/gc_tests.fna");
        let vector = create_normalized_profile(1, FastaNucltudiesIterator::new(dna_test_path.as_path(), 510), &false).1.unwrap();
        assert_eq!(vector[0], 0.6976744186046512);
        assert_eq!(vector[1], 1.0 - 0.6976744186046512);
        assert_eq!(create_printed_vector(1, &ALPHABET_MAP, &REV_ALPHABET_VALUES), vec!["A", "C"])
    }

    #[test]
    fn two_mer() {
        let dna_test_path = PathBuf::from("../tests/gc_tests.fna");

        let vector = create_normalized_profile(2, FastaNucltudiesIterator::new(dna_test_path.as_path(), 512), &false).1.unwrap();
        for (i, v) in [0.2682926829268293, 0.04878048780487805, 0.07317073170731707, 0.17073170731707318, 0.14634146341463414, 0.04878048780487805, 0.024390243902439025, 0.024390243902439025, 0.07317073170731707, 0.12195121951219512].iter().enumerate() {
            assert_eq!(vector[i], v.clone());
        }

        let labels = create_printed_vector(2, &ALPHABET_MAP, &REV_ALPHABET_VALUES);
        for (i, c) in ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"].iter().enumerate() {
            assert_eq!(labels[i], *c);
        }

        let counts = minimize_vector(&mut get_fasta_vector(2, FastaNucltudiesIterator::new(&dna_test_path, 512), &ALPHABET_MAP, &ALPHABET_VALUES).unwrap(), 2);
        for (i, c) in [11.0, 2.0, 3.0, 7.0, 6.0, 2.0, 1.0, 1.0, 3.0, 5.0].iter().enumerate() {
            assert_eq!(counts[i], c.clone());
        }
    }

    #[test]
    fn rev_complement() {
        // k = 1
        for bits in ALPHABET_VALUES.values() {
            let p =  bits.clone();
            let q = bits ^ usize::MAX;
            assert_eq!(p, complement_index(q, 1));
        }

        // k=2
        for bits_a in ALPHABET_VALUES.values() {
            for bits_b in ALPHABET_VALUES.values() {
                let p = (bits_a << 2) | (bits_b);
                let q = ((bits_a) | (bits_b << 2)) ^ usize::MAX;
                assert_eq!(p, complement_index(q, 2));
           }
        }

        // k=3
        for bits_a in ALPHABET_VALUES.values() {
            for bits_b in ALPHABET_VALUES.values() {
                for bits_c in ALPHABET_VALUES.values() {
                    let p = (bits_a << 4) | (bits_b << 2) | bits_c;
                    let q = ((bits_a) | (bits_b << 2) | (bits_c << 4)) ^ usize::MAX;
                    assert_eq!(p, complement_index(q, 3));
                }
            }
        }
    }
}