use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::hash::{Hash, Hasher};
use std::mem::discriminant;
use std::collections::HashMap;
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

impl Display for TaxonomicRank {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxonomicRank::Domain => write!(f, "Domain"),
            TaxonomicRank::Phylum => write!(f, "Phylum"),
            TaxonomicRank::Class => write!(f, "Class"),
            TaxonomicRank::Order => write!(f, "Order"),
            TaxonomicRank::Family => write!(f, "Family"),
            TaxonomicRank::Genus => write!(f, "Genus"),
            TaxonomicRank::Species => write!(f, "Species"),
        }
    }
}
impl Debug for TaxonomicRank {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxonomicRank::Domain => write!(f, "Domain"),
            TaxonomicRank::Phylum => write!(f, "Phylum"),
            TaxonomicRank::Class => write!(f, "Class"),
            TaxonomicRank::Order => write!(f, "Order"),
            TaxonomicRank::Family => write!(f, "Family"),
            TaxonomicRank::Genus => write!(f, "Genus"),
            TaxonomicRank::Species => write!(f, "Species"),
        }
    }
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

    fn taxonomic_rank2prefix(taxonomic_rank: TaxonomicRank) -> char {
        match taxonomic_rank {
            TaxonomicRank::Domain => 'd',
            TaxonomicRank::Phylum => 'p',
            TaxonomicRank::Class => 'c',
            TaxonomicRank::Order => 'o',
            TaxonomicRank::Family => 'f',
            TaxonomicRank::Genus => 'g',
            TaxonomicRank::Species => 's',
        }
    }

    fn get_father(&self) -> Option<TaxonomicRank> {
        match self {
            TaxonomicRank::Domain => None,
            TaxonomicRank::Phylum => Some(TaxonomicRank::Domain),
            TaxonomicRank::Class => Some(TaxonomicRank::Phylum),
            TaxonomicRank::Order => Some(TaxonomicRank::Class),
            TaxonomicRank::Family => Some(TaxonomicRank::Order),
            TaxonomicRank::Genus => Some(TaxonomicRank::Family),
            TaxonomicRank::Species => Some(TaxonomicRank::Genus),
        }
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

#[derive(Serialize, Deserialize, Clone, PartialEq)]
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

    pub fn get_taxa(&self, rank: &TaxonomicRank) -> Option<&String> { self.rank2name.get(rank) }

    pub fn equal_to_rank(&self, other: &Taxonomy, rank: &TaxonomicRank) -> bool {
        match self.lcu(other) {
            None => false,
            Some(r) => r <= *rank,
        }
    }
    
    // Returns the stem of the taxa, up-to, and including, the given rank.
    pub fn limit2rank(&self, rank: &TaxonomicRank) -> Self {
        let rank2name: HashMap<TaxonomicRank, String> = TaxonomicRank::down_iterator()
            .filter_map(|rank_index| if rank_index >= *rank {
                self.rank2name.get(&rank_index).map(|name| (rank_index, name.clone()))
            } else {
                None
            })
            .collect();

        Self {
            rank2name,
        }
    }
    
    // Returns false if the rank is not-named / empty / non-existing.
    pub fn has_rank(&self, rank: &TaxonomicRank) -> bool {
        match self.get_taxa(rank) {
            None => false,
            Some(s) => !s.is_empty(),
        }
    }

    pub fn fill_rank(&mut self, rank: &TaxonomicRank) {
        let father_name = match rank.get_father() {
            Some(father) => self.get_taxa(&father).cloned().unwrap_or(format!("Unknown{}", father)),
            None => "UnknownRank".to_string(),
        };
        
        let new_name = format!("{}->Unknown{}", father_name, rank);
        self.rank2name.insert(*rank, new_name);
    }

}

impl Hash for Taxonomy {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for s in  TaxonomicRank::down_iterator().map(|tr| self.get_taxa(&tr)) {
            if let Some(s) = s {
                s.hash(state);
            } else {
                break
            }
        }
    }
}

impl PartialOrd for Taxonomy {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Taxonomy {
    fn cmp(&self, other: &Self) -> Ordering {
        for tr in TaxonomicRank::down_iterator() {
            let self_taxa = self.get_taxa(&tr);
            let other_taxa = other.get_taxa(&tr);
            if let (Some(st), Some(ot)) = (self_taxa, other_taxa) {
                match st.cmp(ot) {
                    Ordering::Equal => continue,
                    Ordering::Greater => return Ordering::Greater,
                    Ordering::Less => return Ordering::Less,
                }
            } else if self_taxa.is_some() {
                return Ordering::Greater;
            } else if other_taxa.is_some() {
                return Ordering::Less;
            } else {
                return Ordering::Equal;
            }
        }
        Ordering::Equal // Ran through all steps, all were equal.
    }
}

impl Eq for Taxonomy {

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

impl Display for Taxonomy {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        TaxonomicRank::down_iterator()
            .filter_map(|tr| self.get_taxa(&tr).map(|name| if tr != TaxonomicRank::Species {
                format!("{}__{};", TaxonomicRank::taxonomic_rank2prefix(tr), name)
            } else {
                format!("{}__{}", TaxonomicRank::taxonomic_rank2prefix(tr), name)
            }))
            .try_for_each(|tr_name| write!(f, "{}", tr_name))
    }
}

impl Debug for Taxonomy {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.rank2name)
    }
}

#[cfg(test)]
mod tests {
    use crate::taxonomy::{TaxonomicRank, Taxonomy};

    #[test]
    fn test_taxonomy() {
        let raw_string = "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Peptacetobacter;s__Peptacetobacter";
        let taxonomy = Taxonomy::from(raw_string);
        assert_eq!(format!("{}", taxonomy), raw_string);
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Species), Some(&"Peptacetobacter".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Genus), Some(&"Peptacetobacter".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Family), Some(&"Peptostreptococcaceae".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Order), Some(&"Peptostreptococcales".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Class), Some(&"Clostridia".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Phylum), Some(&"Firmicutes_A".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Domain), Some(&"Bacteria".to_string()));
    }

    #[test]
    fn test_taxonomy_space() {
        let raw_string = "d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__MX-02;s__MX-02 sp006954405";
        let taxonomy = Taxonomy::from(raw_string);
        assert_eq!(format!("{}", taxonomy), raw_string);
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Species), Some(&"MX-02 sp006954405".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Genus), Some(&"MX-02".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Family), Some(&"Methanomethylophilaceae".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Order), Some(&"Methanomassiliicoccales".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Class), Some(&"Thermoplasmata".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Phylum), Some(&"Thermoplasmatota".to_string()));
        assert_eq!(taxonomy.get_taxa(&TaxonomicRank::Domain), Some(&"Archaea".to_string()));
    }
}