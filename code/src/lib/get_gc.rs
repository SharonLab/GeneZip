// Returns GC%, a value between 0 and 100
pub fn calc_gc<I>(genome: I) -> f64
where I: IntoIterator<Item=u8> {
    let mut gc_counter = 0_f64;
    let mut length_counter = 0_f64;

    for p in genome {
        match p {
            c if [b'A', b'T', b'G', b'C'].contains(&c) => {
                length_counter += 1.0;
                if [b'G', b'C'].contains(&c) {
                    gc_counter += 1.0;
                }
            },
            _ => {},
        }
    }

    if length_counter > 0_f64 {
        (gc_counter / length_counter) * 100.0
    } else {
        0_f64
    }
}


#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;
    use crate::get_gc::calc_gc;

    #[test]
    fn gc_test_file() {
        let dna_test_path = PathBuf::from("../tests/gc_tests.fna");
        assert_eq!(calc_gc(FastaNucltudiesIterator::new(&dna_test_path, 10)), 30.23255813953488);
    }
}