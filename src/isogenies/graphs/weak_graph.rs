use std::collections::{HashSet, VecDeque};
use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId};

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    pub(crate) fn weak_neighbors(&self, node: IsogenyGraphNodeId) -> Vec<IsogenyGraphNodeId> {
        let mut neighbors = HashSet::new();
        neighbors.extend(self.outgoing_edges(node).map(|edge| edge.target()));
        neighbors.extend(self.incoming_edges(node).map(|edge| edge.source()));

        let mut neighbors = neighbors.into_iter().collect::<Vec<_>>();
        neighbors.sort();
        neighbors
    }

    pub(crate) fn weak_bfs_distances(&self, root: IsogenyGraphNodeId) -> Vec<Option<usize>> {
        let mut distances = vec![None; self.node_count()];
        let mut queue = VecDeque::from([root]);
        distances[root.0] = Some(0);

        while let Some(current) = queue.pop_front() {
            let current_distance =
                distances[current.0].expect("queued weak-BFS nodes should already have a distance");
            for neighbor in self.weak_neighbors(current) {
                if distances[neighbor.0].is_none() {
                    distances[neighbor.0] = Some(current_distance + 1);
                    queue.push_back(neighbor);
                }
            }
        }

        distances
    }

    pub(crate) fn weak_degree(&self, node: IsogenyGraphNodeId) -> usize {
        self.weak_neighbors(node).len()
    }

    pub(crate) fn has_self_loop(&self, node: IsogenyGraphNodeId) -> bool {
        self.outgoing_edges(node).any(|edge| edge.target() == node)
    }
}
