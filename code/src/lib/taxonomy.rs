use std::cmp::Ordering;
use std::mem::discriminant;
use hashbrown::HashMap;
use serde::{Serialize, Deserialize};


#[derive(Serialize, Deserialize, Clone, Copy, Hash, PartialEq)]
pub enum TaxonomicRank {
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
}

impl TaxonomicRank {
    #[allow(dead_code)]
    fn is_same_rank(&self, other: &TaxonomicRank) -> bool {
        discriminant(self) == discriminant(other)
    }

    fn rank2distance(&self) -> usize {
        match self {
            TaxonomicRank::Domain => 6,
            TaxonomicRank::Phylum => 5,
            TaxonomicRank::Class => 4,
            TaxonomicRank::Order => 3,
            TaxonomicRank::Family => 2,
            TaxonomicRank::Genus => 1,
            TaxonomicRank::Species => 0,
        }
    }

    fn down_iterator() -> impl Iterator<Item = TaxonomicRank> {
        [TaxonomicRank::Domain, TaxonomicRank::Phylum, TaxonomicRank::Class, TaxonomicRank::Order, TaxonomicRank::Family, TaxonomicRank::Genus, TaxonomicRank::Species].iter().copied()
    }
}

impl Ord for TaxonomicRank {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rank2distance().cmp(&other.rank2distance())
    }
}

impl PartialOrd for TaxonomicRank {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }

    fn lt(&self, other: &Self) -> bool {
        self.rank2distance().lt(&other.rank2distance())
    }

    fn le(&self, other: &Self) -> bool {
        self.rank2distance().le(&other.rank2distance())
    }

    fn gt(&self, other: &Self) -> bool {
        self.rank2distance().gt(&other.rank2distance())
    }

    fn ge(&self, other: &Self) -> bool {
        self.rank2distance().ge(&other.rank2distance())
    }
}

impl Eq for TaxonomicRank {

}

#[derive(Serialize, Deserialize, Clone)]
pub struct Taxonomy {
    rank2name: HashMap<TaxonomicRank, String>
}

impl Taxonomy {
    pub fn lcu(&self, other: &Taxonomy) -> Option<TaxonomicRank> {
        let mut found_lcu = None;
        for rank in TaxonomicRank::down_iterator() {
            if let (Some(ns), Some(no)) = (self.rank2name.get(&rank), other.rank2name.get(&rank)) {
                if ns == no {
                    found_lcu = Some(rank);
                } else {
                    break
                }
            } else {
                break
            }
        }

        found_lcu
    }

    #[allow(dead_code)]
    pub fn get_taxa(&self, rank: &TaxonomicRank) -> Option<&String> { self.rank2name.get(rank) }

    pub fn equal_to_rank(&self, other: &Taxonomy, rank: &TaxonomicRank) -> bool {
        match self.lcu(other) {
            None => false,
            Some(r) => r <= *rank,
        }
    }
}

fn taxa_parser(value: &str) -> (TaxonomicRank, &str) {
    if value.len() < 3 {
        panic!("E: invalid taxonomic rank! got '{}', quitting", value);
    }
    let name = &value[3..];
    match value.as_bytes()[0] {
        b'd' => (TaxonomicRank::Domain, name),
        b'p' => (TaxonomicRank::Phylum, name),
        b'c' => (TaxonomicRank::Class, name),
        b'o' => (TaxonomicRank::Order, name),
        b'f' => (TaxonomicRank::Family, name),
        b'g' => (TaxonomicRank::Genus, name),
        b's' => (TaxonomicRank::Species, name),
        _ => panic!("E: invalid taxonomy, got '{}', it has no rank", value),
    }
}

impl From<&str> for Taxonomy {
    fn from(value: &str) -> Self {
        let mut rank2name = HashMap::new();

        for rank in value.split(';') {
            let (rank, name) = taxa_parser(rank);
            rank2name.insert(rank,name.to_string());
        }

        Taxonomy {
            rank2name
        }
    }
}