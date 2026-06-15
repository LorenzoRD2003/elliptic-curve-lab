use std::collections::VecDeque;
use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId};

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Returns the weakly connected components of the current directed graph.
    ///
    /// This treats every stored edge `u -> v` as an undirected adjacency for the
    /// purpose of component discovery. Components and the node ids inside each
    /// component are returned in a stable dense-id discovery order.
    pub fn weakly_connected_components(&self) -> Vec<Vec<IsogenyGraphNodeId>> {
        let mut visited = vec![false; self.node_count()];
        let mut components = Vec::new();

        for node in self.nodes() {
            let start = node.id();
            if visited[start.0] {
                continue;
            }

            let mut component = Vec::new();
            let mut queue = VecDeque::from([start]);
            visited[start.0] = true;

            while let Some(current) = queue.pop_front() {
                component.push(current);
                for neighbor in self.weak_neighbors(current) {
                    if !visited[neighbor.0] {
                        visited[neighbor.0] = true;
                        queue.push_back(neighbor);
                    }
                }
            }
            components.push(component);
        }
        components
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Fp, traits::Field};
    use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};

    type F41 = Fp<41>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn weak_components_of_depth_zero_graph_are_one_isolated_node() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        assert_eq!(
            graph.weakly_connected_components(),
            vec![vec![IsogenyGraphNodeId(0)]]
        );
    }

    #[test]
    fn weak_components_ignore_edge_direction_for_single_edge_graph() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        assert_eq!(
            graph.weakly_connected_components(),
            vec![vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]]
        );
    }

    #[test]
    fn weak_components_stay_stable_for_larger_connected_example() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(2)
            .build()
            .expect("depth-two graph should build");

        let components = graph.weakly_connected_components();

        assert_eq!(components.len(), 1);
        assert_eq!(components[0][0], IsogenyGraphNodeId(0));
        assert_eq!(components[0].len(), graph.node_count());
    }
}
