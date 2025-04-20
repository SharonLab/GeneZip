use super::graph::Graph;
use std::collections::{HashSet, HashMap};
use std::cmp::Ordering;
use serde::{Serialize, Deserialize};

#[derive(Clone, Serialize, Deserialize)]
struct Edge<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> {
    #[serde(bound(serialize = "VertexType: Serialize,\
    EdgeData: Serialize",
    deserialize = "VertexType: Deserialize<'de>,\
                  EdgeData: Deserialize<'de>"))]
    other_vertex: VertexType,
    data: EdgeData, // https://users.rust-lang.org/t/how-to-serialize-deserialize-an-async-std-rwlock-t-where-t-serialize-deserialize/37407/2
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> Edge<VertexType, EdgeData> {
    fn new(other_vertex: &VertexType, data: EdgeData) -> Self {
        Edge {
            other_vertex: (*other_vertex).clone(),
            data,
        }
    }
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> PartialOrd for Edge<VertexType, EdgeData> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }

    fn lt(&self, other: &Self) -> bool {
        self.other_vertex < other.other_vertex
    }

    fn le(&self, other: &Self) -> bool {
        self.other_vertex <= other.other_vertex
    }

    fn gt(&self, other: &Self) -> bool {
        self.other_vertex > other.other_vertex
    }

    fn ge(&self, other: &Self) -> bool {
        self.other_vertex >= other.other_vertex
    }
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> PartialEq for Edge<VertexType, EdgeData> {
    fn eq(&self, other: &Self) -> bool {
        self.other_vertex == other.other_vertex
    }
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> Eq for Edge<VertexType, EdgeData> {
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> Ord for Edge<VertexType, EdgeData> {
    fn cmp(&self, other: &Self) -> Ordering {
        if self < other {
            Ordering::Less
        } else if self == other {
            Ordering::Equal
        } else {
            Ordering::Greater
        }
    }

    fn max(self, other: Self) -> Self where Self: Sized {
        if self > other {
            self
        } else {
            other
        }
    }

    fn min(self, other: Self) -> Self where Self: Sized {
        if self < other {
            self
        } else {
            other
        }
    }
}


#[derive(Clone, Serialize, Deserialize)]
struct VertexEdges<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> {
    base_vertex: VertexType,
    edges: Vec<Edge<VertexType, EdgeData>>,  // Vec<(VertexType, Option<EdgeData>)>
}

impl<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone> VertexEdges<VertexType, EdgeData> {
    fn new(base_vertex: VertexType) -> Self {
        VertexEdges {
            base_vertex,
            edges: Vec::new(),
        }
    }
}

#[derive(Clone, Serialize, Deserialize)] // debug
pub struct StaticGraph<VertexType: std::clone::Clone + std::cmp::PartialOrd + std::cmp::Eq, EdgeData: std::clone::Clone, GraphMetadata> {
    edge_data: Vec<Option<EdgeData>>,
    edges: Vec<VertexEdges<VertexType, usize>>,
    metadata: GraphMetadata,
    number_of_edges: usize,
}

impl<VertexType, EdgeData, GraphMetadata> StaticGraph<VertexType, EdgeData, GraphMetadata>
    where VertexType: std::cmp::Eq,
          VertexType: std::hash::Hash,
          VertexType: std::clone::Clone,
          VertexType: std::fmt::Display,
          VertexType: std::cmp::PartialOrd,
          VertexType: std::cmp::Ord,
          VertexType: std::marker::Send,
          VertexType: std::marker::Sync,
          EdgeData: std::clone::Clone,
          EdgeData: std::marker::Send,
          EdgeData: std::marker::Sync,
          GraphMetadata: std::clone::Clone,
          GraphMetadata: std::marker::Send,
          GraphMetadata: std::marker::Sync, {

    // The code is wasting too much time in get_children which is called a lot
    // A solution is to keep the edge mention on both vertices.
    // To do so, means we need to create the same edge twice, once a -> b and once b -> a
    // This helper does a -> b, and can be called twice.
    //
    // This code should only happen if the edges are being added, not if the data is being updated.
    fn add_edge_helper(&mut self, edge_from: &VertexType, edge_to: &VertexType, data: usize) -> bool {
        match self.edges.binary_search_by_key(&edge_from, |v| &v.base_vertex) {
            Ok(edge_from_index) => {
                // On failure to merge with existing data, create a new edge
                self.edges[edge_from_index].edges.push(Edge::new(edge_to, data));
                self.edges[edge_from_index].edges.sort();
                true
            },
            Err(_) => false,
        }
    }

    // Same as add_edge_helper, but for dropping the edge
    fn drop_edge_helper(&mut self, edge_from: &VertexType, edge_to: &VertexType) -> Option<EdgeData> {
        match self.edges.binary_search_by_key(&edge_from, |v| &v.base_vertex) {
            Ok(edge_from_index) => match self.edges[edge_from_index].edges.binary_search_by_key(&edge_to, |v| &v.other_vertex) {
                Ok(edge_to_index) => {
                    let data_index = self.edges[edge_from_index].edges.remove(edge_to_index).data;
                    self.edge_data[data_index].take()
                },
                _ => None,
            },
            _ => None,
        }
    }
}

impl<VertexType, EdgeData, GraphMetadata> Graph for StaticGraph<VertexType, EdgeData, GraphMetadata>
    where VertexType: std::cmp::Eq,
          VertexType: std::hash::Hash,
          VertexType: std::clone::Clone,
          VertexType: std::fmt::Display,
          VertexType: std::cmp::PartialOrd,
          VertexType: std::cmp::Ord,
          VertexType: std::marker::Send,
          VertexType: std::marker::Sync,
          EdgeData: std::clone::Clone,
          EdgeData: std::marker::Send,
          EdgeData: std::marker::Sync,
          GraphMetadata: std::clone::Clone,
          GraphMetadata: std::marker::Send,
          GraphMetadata: std::marker::Sync, {
    type VertexType = VertexType;
    type EdgeData = EdgeData;
    type GraphMetadata = GraphMetadata;

    fn new(number_of_expected_vertices: Option<usize>, metadata: GraphMetadata) -> Self {
        StaticGraph {
            edge_data: Vec::new(),
            edges: match number_of_expected_vertices {
                Some(n) => Vec::with_capacity(n),
                None => Vec::new(),
            },
            metadata,
            number_of_edges: 0,
        }
    }

    // Can only be used on an empty graph
    // Duplicate vertices will be removed
    fn batch_add_vertex(&mut self, vertices: &[VertexType], edges: &[(VertexType, VertexType, Option<EdgeData>)], edge_collision: fn(Option<Self::EdgeData>, Option<Self::EdgeData>) -> Option<Self::EdgeData>) {
        if self.edges.is_empty() {
            for vertex in vertices {
                self.edges.push(VertexEdges::new(vertex.clone()));
            }
            self.edges.sort_by(|u, v| u.base_vertex.cmp(&v.base_vertex));
            self.edges.dedup_by(|u, v| u.base_vertex.eq(&v.base_vertex));

            for (u, v, data) in edges {
                self.add_edge(u, v, &data, edge_collision);
            }
        }
    }

    fn mut_metadata(&mut self) -> &mut Self::GraphMetadata {
        &mut self.metadata
    }

    fn metadata(&self) -> &Self::GraphMetadata {
        &self.metadata
    }

    fn get_vertex(&self, vertex: &Self::VertexType) -> Option<&Self::VertexType> {
        match self.edges.binary_search_by_key(&vertex, |v| &v.base_vertex) {
            Ok(index) => Some(&self.edges[index].base_vertex),
            Err(_) => None,
        }
    }

    // fn has_vertex(&self, vertex: &Self::VertexType) -> bool {
    //     self.edges.binary_search_by_key(&vertex, |v| &v.base_vertex).is_ok()
    // }

    fn add_vertex(&mut self, vertex: &Self::VertexType) -> bool {
        if self.has_vertex(vertex) {
            true
        } else {
            self.edges.push(VertexEdges::new(vertex.clone()));
            self.edges.sort_unstable_by(|u, v| u.base_vertex.cmp(&v.base_vertex));
            true
        }
    }

    fn add_edge(&mut self, edge_from: &Self::VertexType, edge_to: &Self::VertexType, data: &Option<Self::EdgeData>, edge_collision: fn(Option<Self::EdgeData>, Option<Self::EdgeData>) -> Option<Self::EdgeData>) -> bool {
        if !(self.has_vertex(&edge_from) && self.has_vertex(&edge_to)) {
            false
        } else {
            match self.edges.binary_search_by_key(&edge_from, |v| &v.base_vertex) {
                Ok(edge_from_index) => {
                    // Try to merge with existing data
                    if let Ok(edge_to_index) = self.edges[edge_from_index].edges.binary_search_by_key(&edge_to, |v| &v.other_vertex) {
                        let existing_data_index = self.edges[edge_from_index].edges[edge_to_index].data;
                        let merged_data = (edge_collision)(self.edge_data.remove(existing_data_index), data.clone());
                        self.edge_data.insert(existing_data_index, merged_data);
                    } else {
                        // On failure to merge with existing data, create a new edge
                        self.edge_data.push(data.clone());
                        let data_index = self.edge_data.len() - 1;
                        if edge_from != edge_to {
                            if self.add_edge_helper(edge_from, edge_to, data_index) !=
                                self.add_edge_helper(edge_to, edge_from, data_index) {
                                panic!("ERROR: static_graph: add_edge: out of balance, can't recover");
                            }
                        } else {
                            self.add_edge_helper(edge_from, edge_to, data_index);
                        }
                    }
                    self.number_of_edges += 1;
                    true
                },
                Err(_) => false,
            }
        }
    }

    fn has_edge(&self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> bool {
        if let Ok(edge_from_index) = self.edges.binary_search_by_key(&edge_from, |v| &v.base_vertex) {
            self.edges[edge_from_index].edges.binary_search_by_key(&edge_to, |v| &v.other_vertex).is_ok()
        } else {
            false
        }
    }

    fn get_edge_data(&self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> &Option<Self::EdgeData> {
        match self.edges.binary_search_by_key(&edge_from, |v| &v.base_vertex) {
            Ok(edge_from_index) => match self.edges[edge_from_index].edges.binary_search_by_key(&edge_to, |v| &v.other_vertex) {
                Ok(edge_to_index) => &self.edge_data[self.edges[edge_from_index].edges[edge_to_index].data],
                _ => &None,
            },
            _ => &None,
        }
    }

    fn get_edges_from(&self, vertex: &Self::VertexType) -> HashMap<&Self::VertexType, &Option<Self::EdgeData>> {
        let mut result = HashMap::new();
        if let Ok(vertex_edges_index) = self.edges.binary_search_by_key(&vertex, |v| &v.base_vertex) {
            for (other_vertex, data) in self.edges[vertex_edges_index].edges.iter().map(|an_edge| (&an_edge.other_vertex, &an_edge.data)) {
               result.insert(other_vertex, &self.edge_data[*data]);
            }
        }
        result
    }

    fn get_children(&self, vertex: &Self::VertexType) -> HashSet<&Self::VertexType> {
        let mut result = HashSet::new();
        if let Ok(vertex_edges_index) = self.edges.binary_search_by_key(&vertex, |v| &v.base_vertex) {
            for edge in &self.edges[vertex_edges_index].edges {
                result.insert(&edge.other_vertex);
            }
        }
        result
    }

    fn drop_edge(&mut self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> Option<Self::EdgeData> {
        if self.has_edge(edge_from, edge_to) {
            self.number_of_edges -= 1;
            let results = self.drop_edge_helper(edge_to, edge_from);
            if edge_from != edge_to {
                self.drop_edge_helper(edge_from, edge_to);
            }
            results
        } else {
            None
        }
    }

    fn drop_vertex(&mut self, vertex: &Self::VertexType) {
        if let Ok(vertex_index) = self.edges.binary_search_by_key(&vertex, |v| &v.base_vertex) {
            for edge_to_other_vertex in &self.edges.remove(vertex_index).edges { // Remove us, then removes us from all the vertices we have edges with
                self.number_of_edges -= 1;
                if vertex != &edge_to_other_vertex.other_vertex {
                    if let Ok(other_vertex_index) = self.edges.binary_search_by_key(&&edge_to_other_vertex.other_vertex, |v| &v.base_vertex) {
                        let other_vertex = &mut self.edges[other_vertex_index];
                        let mut to_remove = None;
                        // Find the other vertex edge with us
                        for (loop_index, other_vertex_edge_data) in other_vertex.edges.iter().enumerate() {
                            if &other_vertex_edge_data.other_vertex == vertex {
                                to_remove = Some(loop_index);
                                break;
                            }
                        }
                        // Remove it
                        if let Some(to_remove_index) = to_remove {
                            other_vertex.edges.remove(to_remove_index);
                        }
                    }
                }
            }
        }
    }

    fn vertices(&self) -> Vec<&Self::VertexType> {
        let mut results = Vec::with_capacity(self.edges.len());
        self.edges.iter().for_each(|vertex| results.push(&vertex.base_vertex));
        results
    }

    fn len(&self) -> usize { self.edges.len() }

    fn edge_triplets(&self) -> Vec<(&Self::VertexType, &Self::VertexType, &Option<Self::EdgeData>)> {
        let mut results: Vec<(&Self::VertexType, &Self::VertexType, &Option<Self::EdgeData>)> = Vec::with_capacity(self.number_of_edges());
        let mut temp_holder = Vec::with_capacity(self.number_of_edges());
        for edges_from in self.edges.iter() {
            for (edge_to, data) in edges_from.edges.iter().map(|an_edge| (&an_edge.other_vertex, &an_edge.data)) {
                if &edges_from.base_vertex < edge_to { // Each edge is held twice, but we want to add it only once.
                    temp_holder.push((&edges_from.base_vertex, edge_to, &self.edge_data[*data]));
                }
            }
        }

        for (edge_from, edge_to, v) in temp_holder.iter() {
            results.push((*edge_from, *edge_to, v));
        }

        results
    }

    fn number_of_edges(&self) -> usize {
        self.number_of_edges
    }
}

#[cfg(test)]
mod tests {
    use crate::static_graph::StaticGraph;
    use crate::graph::Graph;

    #[test]
    fn graph_structure() {
        let mut graph: StaticGraph<u8, u8, Option<u8>> = StaticGraph::new(Some(10), None);

        assert_eq!(graph.add_vertex(&1), true);
        assert_eq!(graph.add_vertex(&2), true);
        assert_eq!(graph.add_vertex(&3), true);
        assert_eq!(graph.add_vertex(&1), true);
        assert_eq!(graph.add_vertex(&2), true);
        assert_eq!(graph.add_vertex(&3), true);

        assert_eq!(graph.has_vertex(&1), true);
        assert_eq!(graph.has_vertex(&2), true);
        assert_eq!(graph.has_vertex(&3), true);
        assert_eq!(graph.has_vertex(&4), false);

        assert_eq!(graph.add_edge(&1, &1, &Some(2), |a, _| a), true);
        assert_eq!(graph.add_edge(&1, &2, &Some(3), |a, _| a), true);
        assert_eq!(graph.add_edge(&3, &2, &Some(5), |a, _| a), true);
        assert_eq!(graph.add_edge(&4, &2, &Some(6), |a, _| a), false);


        assert_eq!(graph.has_edge(&1, &1), true);
        assert_eq!(graph.has_edge(&1, &2), true);
        assert_eq!(graph.has_edge(&2, &2), false);
        assert_eq!(graph.has_edge(&1, &3), false);
        assert_eq!(graph.has_edge(&3, &2), true);
        assert_eq!(graph.has_edge(&4, &2), false);

        assert_eq!(graph.len(), 3);

        assert_eq!(graph.number_of_edges(), 3);

        // Test get children
        let children_of_one = graph.get_children(&1);
        assert_eq!(children_of_one.len(), 2);
        assert_eq!(children_of_one.contains(&1), true);
        assert_eq!(children_of_one.contains(&2), true);

        // Test get_connected_component
        for num in 1..=3 {
            let connected_component = graph.get_connected_component(&num);
            for i in 1..=3 {
                println!("static_graph: graph_structure: {} in graph.get_connected_component(&{})", i, num);
                assert_eq!(connected_component.contains(&i), true);
            }
        }

        // connectivity_component_representatives
        let connectivity_component_representatives = graph.get_connected_component_representatives();
        assert_eq!(connectivity_component_representatives.contains(&1), true);
        assert_eq!(connectivity_component_representatives.contains(&2), false);
        assert_eq!(connectivity_component_representatives.contains(&3), false);

        assert_eq!(graph.drop_edge(&1, &1), Some(2));
        assert_eq!(graph.drop_edge(&1, &2), Some(3));
        assert_eq!(graph.number_of_edges(), 1);
        assert_eq!(graph.has_edge(&1, &1), false);
        assert_eq!(graph.has_edge(&1, &2), false);

        graph.drop_vertex(&1);
        assert_eq!(graph.len(), 2);
        assert_eq!(graph.number_of_edges(), 1);
        assert_eq!(graph.has_vertex(&1), false);

        graph.drop_vertex(&3);
        assert_eq!(graph.len(), 1);
        assert_eq!(graph.number_of_edges(), 0);
        assert_eq!(graph.has_edge(&3, &2), false);
        assert_eq!(graph.has_vertex(&3), false);
    }
}