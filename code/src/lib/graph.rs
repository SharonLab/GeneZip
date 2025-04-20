use std::collections::{HashMap, HashSet, LinkedList};

fn bfs_backtrack_tree2path<'a, VertexType>(backtrack_tree: &HashMap<&VertexType, Option<&'a VertexType>>, node: &'a VertexType) -> Vec<&'a VertexType>
    where VertexType: std::clone::Clone,
          VertexType: std::cmp::Eq,
          VertexType: std::hash::Hash, {
    let mut path: Vec<&VertexType> = vec![node];
    while let Some(Some(backtrack_node)) = backtrack_tree.get(&path[0]) {
        path.insert(0, *backtrack_node);
    }
    path
}

pub trait Graph: Send {
    type VertexType: std::cmp::Eq + std::hash::Hash + std::clone::Clone + std::fmt::Display + std::cmp::PartialOrd + std::marker::Send;
    type EdgeData: std::clone::Clone + std::marker::Send;
    type GraphMetadata: std::clone::Clone + std::marker::Send;

    fn new(number_of_expected_vertices: Option<usize>, metadata: Self::GraphMetadata) -> Self;

    fn batch_add_vertex(&mut self, vertices: &[Self::VertexType], edges: &[(Self::VertexType, Self::VertexType, Option<Self::EdgeData>)], edge_collision: fn(Option<Self::EdgeData>, Option<Self::EdgeData>) -> Option<Self::EdgeData>);

    fn mut_metadata(&mut self) -> &mut Self::GraphMetadata;
    fn metadata(&self) -> &Self::GraphMetadata;

    // Returns the internal reference to an equal vertex.
    fn get_vertex(&self, vertex: &Self::VertexType) -> Option<&Self::VertexType>;

    fn has_vertex(&self, vertex: &Self::VertexType) -> bool {
        self.get_vertex(vertex).is_some()
    }

    // |V|
    fn len(&self) -> usize;

    //
    // Returns true if the vertex now exists in the graph, not if the vertex was added.
    // Adding the same vertex twice will return true both times.
    //
    fn add_vertex(&mut self, vertex: &Self::VertexType) -> bool;

    // Returns true if the edge exists (added, or pre existing), other wise false
    fn add_edge(&mut self, edge_from: &Self::VertexType, edge_to: &Self::VertexType, data: &Option<Self::EdgeData>, edge_collision: fn(Option<Self::EdgeData>, Option<Self::EdgeData>) -> Option<Self::EdgeData>) -> bool;

    fn has_edge(&self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> bool;

    fn get_edge_data(&self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> &Option<Self::EdgeData>;

    fn get_edges_from(&self, vertex: &Self::VertexType) -> HashMap<&Self::VertexType, &Option<Self::EdgeData>>;

    fn get_children(&self, vertex: &Self::VertexType) -> HashSet<&Self::VertexType>;

    fn drop_edge(&mut self, edge_from: &Self::VertexType, edge_to: &Self::VertexType) -> Option<Self::EdgeData>;

    fn drop_vertex(&mut self, vertex: &Self::VertexType);

    fn vertices(&self) -> Vec<&Self::VertexType>;

    fn edge_triplets(&self) -> Vec<(&Self::VertexType, &Self::VertexType, &Option<Self::EdgeData>)>;

    fn number_of_edges(&self) -> usize;

    // Graph tools:
    fn bfs<'a>(&'a self, start: &'a Self::VertexType, goal: &'a Self::VertexType) -> Option<Vec<&'a Self::VertexType>> {
        let mut queue: Vec<&Self::VertexType> = vec![start];
        let mut visited = HashSet::new();
        visited.insert(start);
        let mut backtrack_tree: HashMap<&Self::VertexType, Option<&Self::VertexType>> = HashMap::new();
        backtrack_tree.insert(start, None);
        while !queue.is_empty() {
            if let Some(node) = queue.pop() {
                if node == goal {
                    return Some(bfs_backtrack_tree2path(&backtrack_tree, &node));
                }
                for child in self.get_children(&node) {
                    if !visited.contains(child) {
                        backtrack_tree.insert(child, Some(node));
                        visited.insert(child);
                        queue.push(child);
                    }
                }
            }
        }
        None
    }

    // Returns the unweighted distance between two vertices.
    // If no path exists, returns None.
    fn unweighted_distance(&self, from: &Self::VertexType, to: &Self::VertexType) -> Option<usize> {
        let path = self.bfs(from, to);
        if let Some(path) = path {
            Some(path.len() - 1)
        } else {
            None
        }
    }

    //
    // Using BFS, checks if the starting point is part of a circle.
    //
    // Returns None or the full circle
    //
    fn check_circle<'a>(&'a self, starting_point: &'a Self::VertexType) -> Option<Vec<&'a Self::VertexType>> {
        let mut path = Vec::new();

        // Edge case of one vertex with self edge
        {
            let children = self.get_children(starting_point);
            if children.len() == 1 && children.contains(starting_point) { // Self edge
                path.push(starting_point);
                return Some(path);
            }
        }

        // Handle the general case
        let mut current_step = starting_point;
        let mut prev_step = starting_point;
        loop {
            path.push(current_step);
            let children = self.get_children(current_step);

            if children.len() != 2 || children.contains(current_step) { // Either a split or end of path OR self edge
                return None;
            }

            let next_step = match children.iter().find(|&vertex| vertex != &prev_step) {
                Some(vertex) => *vertex,
                None => return None, // This should NEVER happen
            };

            if next_step == starting_point { // we finished
                return Some(path);
            }

            prev_step = current_step;
            current_step = next_step;
        }
    }

    fn get_connected_component_representatives(&self) -> HashSet<&Self::VertexType> {
        let mut representatives = HashSet::new(); // A set of vertices that represent connectivity component, one representative per component
        let mut represented: HashSet<Self::VertexType> = HashSet::new(); // A set of all vertices that already show up in one connectivity component

        for vertex in self.vertices() {
            if !represented.contains(&vertex) {
                for part in self.get_connected_component(&vertex).drain() {
                    represented.insert(part.clone());
                }

                representatives.insert(vertex);
            }
        }

        representatives
    }

    fn get_connected_component(&self, vertex: &Self::VertexType) -> HashSet<&Self::VertexType> {
        let mut future_vertices: LinkedList<&Self::VertexType> = LinkedList::new(); // Vertices to visit
        let vertex = self.get_vertex(vertex);
        future_vertices.push_front(vertex.unwrap());
        let mut already_found = HashSet::new(); // Vertices that were visited
        while !future_vertices.is_empty() { // Loop while we find new vertices
            if let Some(vertex) = future_vertices.pop_front() {
                let already_found_len = already_found.len();
                already_found.insert(vertex);
                if already_found.len() == already_found_len + 1 {
                    for child in self.get_children(&vertex) {
                        if !already_found.contains(&child) {
                            future_vertices.push_front(child);
                        }
                    }
                }
            }
        }

        already_found
    }

    fn connected_components(&self) -> Vec<HashSet<&Self::VertexType>> {
        self.get_connected_component_representatives()
            .iter()
            .map(|&r| self.get_connected_component(r))
            .collect()
    }

    fn is_empty(&self) -> bool { self.len() == 0 }
}