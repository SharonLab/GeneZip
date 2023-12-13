pub fn get_contig_name(full_id: &[u8]) -> &[u8] {
    let mut parts = full_id.rsplitn(2, |b| *b == b'_');
    parts.next();
    parts.next().unwrap_or_else(|| panic!("E: Tried to extract contig name from '{:?}', by splitting right-most underline _, however no two parts were found",
                                          full_id))
}

pub fn are_genes_of_same_contig(id1: &[u8], id2: &[u8]) -> bool {
    let id1 = get_contig_name(id1);
    let id2 = get_contig_name(id2);
    id1.eq(id2)
}

pub fn sequence_id2str(id: &[u8]) -> &str {
    std::str::from_utf8(id).unwrap_or_else(|_| panic!("E: Sequence id is not a valid string! this should never happen. The id: '{:?}'", id))
}


#[cfg(test)]
mod tests {
    use crate::contig_naming::{are_genes_of_same_contig, get_contig_name};

    #[test]
    fn test_get_contig_name() {
        let full_name = "test_0_1".as_bytes();
        assert_eq!(get_contig_name(full_name), "test_0".as_bytes());
    }

    #[test]
    fn test_are_genes_of_same_contig() {
        let id_1 = "test_0_1".as_bytes();
        let id_2 = "test_0_2".as_bytes();
        let id_3 = "test_1_2".as_bytes();
        assert!(are_genes_of_same_contig(id_1, id_2));
        assert!(!are_genes_of_same_contig(id_1, id_3));
    }

}