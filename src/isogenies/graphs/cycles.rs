use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId};

/// Returns whether the directed graph contains any directed cycle.
///
/// This includes self-loops as cycles of length `1`.
pub fn has_directed_cycle<C: GraphCurveModel>(graph: &IsogenyGraph<C>) -> bool
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let mut state = vec![VisitState::Unvisited; graph.node_count()];
    for node in graph.nodes() {
        if state[node.id().0] == VisitState::Unvisited
            && dfs_has_cycle(graph, node.id(), &mut state)
        {
            return true;
        }
    }
    false
}

/// Finds simple directed cycles up to `max_len` and returns them in a stable,
/// deduplicated form.
///
/// Each returned cycle is represented as a list of node ids without repeating
/// the start node at the end. Cycles are deduplicated up to cyclic rotation.
/// For example, `v0 -> v1 -> v2 -> v0` is returned once as
/// `[v0, v1, v2]`, not separately from `[v1, v2, v0]`.
pub fn find_small_directed_cycles<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    max_len: usize,
) -> Vec<Vec<IsogenyGraphNodeId>>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    if max_len == 0 {
        return Vec::new();
    }

    let mut search = CycleSearch {
        graph,
        max_len,
        seen: HashSet::new(),
        cycles: Vec::new(),
    };

    for node in graph.nodes() {
        let start = node.id();
        let mut path = vec![start];
        let mut visited = HashSet::from([start]);
        search.dfs_collect(start, start, &mut path, &mut visited);
    }

    search.cycles.sort();
    search.cycles
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum VisitState {
    Unvisited,
    Visiting,
    Visited,
}

fn dfs_has_cycle<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    current: IsogenyGraphNodeId,
    state: &mut [VisitState],
) -> bool
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    state[current.0] = VisitState::Visiting;
    for neighbor in outgoing_neighbors(graph, current) {
        match state[neighbor.0] {
            VisitState::Unvisited => {
                if dfs_has_cycle(graph, neighbor, state) {
                    return true;
                }
            }
            VisitState::Visiting => return true,
            VisitState::Visited => {}
        }
    }
    state[current.0] = VisitState::Visited;
    false
}

fn outgoing_neighbors<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    node: IsogenyGraphNodeId,
) -> Vec<IsogenyGraphNodeId>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let mut neighbors = graph
        .outgoing_edges(node)
        .map(|edge| edge.target())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    neighbors.sort();
    neighbors
}

fn canonical_cycle(cycle: &[IsogenyGraphNodeId]) -> Vec<IsogenyGraphNodeId> {
    if cycle.len() <= 1 {
        return cycle.to_vec();
    }

    let mut best = rotate_cycle(cycle, 0);
    for shift in 1..cycle.len() {
        let candidate = rotate_cycle(cycle, shift);
        if candidate < best {
            best = candidate;
        }
    }
    best
}

fn rotate_cycle(cycle: &[IsogenyGraphNodeId], shift: usize) -> Vec<IsogenyGraphNodeId> {
    cycle[shift..]
        .iter()
        .chain(cycle[..shift].iter())
        .copied()
        .collect()
}

struct CycleSearch<'a, C: GraphCurveModel>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    graph: &'a IsogenyGraph<C>,
    max_len: usize,
    seen: HashSet<Vec<IsogenyGraphNodeId>>,
    cycles: Vec<Vec<IsogenyGraphNodeId>>,
}

impl<C: GraphCurveModel> CycleSearch<'_, C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    fn dfs_collect(
        &mut self,
        start: IsogenyGraphNodeId,
        current: IsogenyGraphNodeId,
        path: &mut Vec<IsogenyGraphNodeId>,
        visited: &mut HashSet<IsogenyGraphNodeId>,
    ) {
        for neighbor in outgoing_neighbors(self.graph, current) {
            if neighbor == start {
                let canonical = canonical_cycle(path);
                if self.seen.insert(canonical.clone()) {
                    self.cycles.push(canonical);
                }
                continue;
            }

            if path.len() >= self.max_len || visited.contains(&neighbor) {
                continue;
            }

            visited.insert(neighbor);
            path.push(neighbor);
            self.dfs_collect(start, neighbor, path, visited);
            path.pop();
            visited.remove(&neighbor);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::{
        IsogenyGraphBuilder, IsogenyGraphNodeId, find_small_directed_cycles, has_directed_cycle,
    };

    type F5 = Fp<5>;
    type F41 = Fp<41>;
    type Curve41 = ShortWeierstrassCurve<F41>;
    type Curve5 = ShortWeierstrassCurve<F5>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn f5_split_two_torsion_curve() -> Curve5 {
        Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
    }

    #[test]
    fn depth_zero_graph_has_no_directed_cycle() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        assert!(!has_directed_cycle(&graph));
        assert!(find_small_directed_cycles(&graph, 4).is_empty());
    }

    #[test]
    fn depth_one_graph_has_no_directed_cycle() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        assert!(!has_directed_cycle(&graph));
        assert!(find_small_directed_cycles(&graph, 4).is_empty());
    }

    #[test]
    fn split_two_torsion_example_has_directed_cycles() {
        let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
            .max_depth(2)
            .build()
            .expect("split two-torsion graph should build");

        assert!(has_directed_cycle(&graph));
    }

    #[test]
    fn self_loop_is_reported_as_length_one_cycle() {
        let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
            .max_depth(2)
            .build()
            .expect("split two-torsion graph should build");

        let cycles = find_small_directed_cycles(&graph, 1);

        assert!(cycles.iter().any(|cycle| cycle.len() == 1));
    }

    #[test]
    fn cycle_search_deduplicates_rotations() {
        let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
            .max_depth(2)
            .build()
            .expect("split two-torsion graph should build");

        let cycles = find_small_directed_cycles(&graph, 3);
        let expected = vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)];

        assert!(cycles.contains(&expected));
        assert_eq!(cycles.iter().filter(|cycle| **cycle == expected).count(), 1);
    }
}
