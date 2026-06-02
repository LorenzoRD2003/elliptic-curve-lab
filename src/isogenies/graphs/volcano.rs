use std::collections::{HashSet, VecDeque};
use std::fmt;
use std::hash::Hash;

use super::{GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId};

/// Educational graph-theoretic role assigned to a node inside one weakly
/// explored component.
///
/// This is intentionally heuristic. It does not attempt to certify the
/// arithmetic meaning of a node, only to attach a small teaching-oriented label
/// derived from weak BFS layering and local weak degree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VolcanoRole {
    /// A node on the first weak-BFS layer that contains reachable
    /// higher-degree branching behavior.
    Surface,
    /// A reachable node with weak degree `1` and no self-loop.
    Floor,
    /// A reachable non-isolated node that is neither classified as
    /// [`Surface`](Self::Surface) nor [`Floor`](Self::Floor).
    Middle,
    /// A reachable node with no weak neighbors at all.
    Isolated,
    /// A node outside the weakly reachable component of the chosen root, or a
    /// node in a graph queried with an invalid root id.
    Unknown,
}

/// Educational weak-BFS layering rooted at one chosen node.
///
/// `levels[i]` contains the reachable nodes at weak distance `i` from the
/// chosen root. Nodes outside that weakly connected component are not placed in
/// any level and are typically labeled [`VolcanoRole::Unknown`].
#[derive(Clone, Debug)]
pub struct VolcanoLayering {
    pub levels: Vec<Vec<IsogenyGraphNodeId>>,
    pub roles: Vec<(IsogenyGraphNodeId, VolcanoRole)>,
}

/// Infers a small volcano-like layering from the weak graph structure around a
/// chosen root.
///
/// This is an educational, graph-theoretic heuristic only. It does not compute
/// endomorphism rings, distinguish ordinary from supersingular components, or
/// prove that the current component is an actual ordinary isogeny volcano in
/// the sense discussed by Sutherland.
///
/// Current heuristic:
///
/// - build weak BFS levels from `root`
/// - compute each reachable node's weak degree, deduplicating reverse edges
/// - mark reachable isolated nodes as [`VolcanoRole::Isolated`]
/// - mark reachable degree-1 nodes without self-loops as
///   [`VolcanoRole::Floor`]
/// - among the remaining reachable higher-degree nodes, find the minimum BFS
///   level and call that shell [`VolcanoRole::Surface`]
/// - mark the other reachable non-isolated nodes as [`VolcanoRole::Middle`]
///
/// Technical debt:
///
/// - The layering depends on the chosen root; it is not an intrinsic crater
///   detection algorithm.
/// - Weak degree is based on distinct weak neighbors and ignores edge
///   multiplicities, so self-loops and parallel edges are only approximated.
/// - The current implementation does not classify edges as ascending,
///   descending, or horizontal.
/// - A true volcanic decomposition would need arithmetic information well
///   beyond this graph-only heuristic.
pub fn infer_volcano_like_layers<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    root: IsogenyGraphNodeId,
) -> VolcanoLayering
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let Some(_) = graph.node(root) else {
        return VolcanoLayering {
            levels: Vec::new(),
            roles: graph
                .nodes()
                .iter()
                .map(|node| (node.id(), VolcanoRole::Unknown))
                .collect(),
        };
    };

    let distances = weak_bfs_distances(graph, root);
    let levels = levels_from_distances(&distances);
    let surface_level = graph
        .nodes()
        .iter()
        .filter_map(|node| {
            let id = node.id();
            distances[id.0].filter(|_| is_high_degree_node(graph, id))
        })
        .min();

    let roles = graph
        .nodes()
        .iter()
        .map(|node| {
            let id = node.id();
            let role = match distances[id.0] {
                None => VolcanoRole::Unknown,
                Some(level) => classify_role(graph, id, level, surface_level),
            };
            (id, role)
        })
        .collect();

    VolcanoLayering { levels, roles }
}

fn classify_role<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    node: IsogenyGraphNodeId,
    level: usize,
    surface_level: Option<usize>,
) -> VolcanoRole
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let total = weak_degree(graph, node);
    if total == 0 {
        VolcanoRole::Isolated
    } else if total == 1 && !has_self_loop(graph, node) {
        VolcanoRole::Floor
    } else if surface_level == Some(level) && is_high_degree_node(graph, node) {
        VolcanoRole::Surface
    } else {
        VolcanoRole::Middle
    }
}

fn weak_bfs_distances<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    root: IsogenyGraphNodeId,
) -> Vec<Option<usize>>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let mut distances = vec![None; graph.node_count()];
    let mut queue = VecDeque::from([root]);
    distances[root.0] = Some(0);

    while let Some(current) = queue.pop_front() {
        let current_distance =
            distances[current.0].expect("queued weak-BFS nodes should already have a distance");
        for neighbor in weak_neighbors(graph, current) {
            if distances[neighbor.0].is_none() {
                distances[neighbor.0] = Some(current_distance + 1);
                queue.push_back(neighbor);
            }
        }
    }

    distances
}

fn levels_from_distances(distances: &[Option<usize>]) -> Vec<Vec<IsogenyGraphNodeId>> {
    let max_level = distances.iter().copied().flatten().max();
    let Some(max_level) = max_level else {
        return Vec::new();
    };

    let mut levels = vec![Vec::new(); max_level + 1];
    for (index, distance) in distances.iter().copied().enumerate() {
        if let Some(level) = distance {
            levels[level].push(IsogenyGraphNodeId(index));
        }
    }
    levels
}

fn weak_degree<C: GraphCurveModel>(graph: &IsogenyGraph<C>, node: IsogenyGraphNodeId) -> usize
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    weak_neighbors(graph, node).len()
}

fn is_high_degree_node<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    node: IsogenyGraphNodeId,
) -> bool
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    weak_degree(graph, node) > 1 || has_self_loop(graph, node)
}

fn has_self_loop<C: GraphCurveModel>(graph: &IsogenyGraph<C>, node: IsogenyGraphNodeId) -> bool
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    graph.outgoing_edges(node).any(|edge| edge.target() == node)
}

fn weak_neighbors<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
    node: IsogenyGraphNodeId,
) -> Vec<IsogenyGraphNodeId>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let mut neighbors = HashSet::new();
    neighbors.extend(graph.outgoing_edges(node).map(|edge| edge.target()));
    neighbors.extend(graph.incoming_edges(node).map(|edge| edge.source()));

    let mut neighbors = neighbors.into_iter().collect::<Vec<_>>();
    neighbors.sort();
    neighbors
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::{
        IsogenyGraph, IsogenyGraphBuilder, IsogenyGraphNodeId, VolcanoRole,
        infer_volcano_like_layers,
    };

    type F17 = Fp<17>;
    type F41 = Fp<41>;
    type Curve17 = ShortWeierstrassCurve<F17>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn graph_with_surface_middle_and_floor_over_f17() -> (IsogenyGraph<Curve17>, IsogenyGraphNodeId)
    {
        for a in 0..17 {
            for b in 0..17 {
                let Ok(curve) = Curve17::new(F17::elem_from_u64(a), F17::elem_from_u64(b)) else {
                    continue;
                };
                let Ok(graph) = IsogenyGraphBuilder::new(curve, 2).max_depth(3).build() else {
                    continue;
                };
                let node_ids = graph
                    .nodes()
                    .iter()
                    .map(|node| node.id())
                    .collect::<Vec<_>>();
                for root in node_ids {
                    let roles = infer_volcano_like_layers(&graph, root)
                        .roles
                        .iter()
                        .map(|(_, role)| *role)
                        .collect::<Vec<_>>();
                    if roles.contains(&VolcanoRole::Surface)
                        && roles.contains(&VolcanoRole::Middle)
                        && roles.contains(&VolcanoRole::Floor)
                    {
                        return (graph, root);
                    }
                }
            }
        }

        panic!("expected to find an F17 example with surface, middle, and floor roles");
    }

    fn role_of(
        layering: &crate::isogenies::graphs::volcano::VolcanoLayering,
        node: IsogenyGraphNodeId,
    ) -> VolcanoRole {
        layering
            .roles
            .iter()
            .find_map(|(candidate, role)| (*candidate == node).then_some(*role))
            .expect("test node should have a role")
    }

    #[test]
    fn invalid_root_marks_all_nodes_unknown() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let layering = infer_volcano_like_layers(&graph, IsogenyGraphNodeId(99));

        assert!(layering.levels.is_empty());
        assert!(
            layering
                .roles
                .iter()
                .all(|(_, role)| *role == VolcanoRole::Unknown)
        );
    }

    #[test]
    fn isolated_root_is_classified_as_isolated() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let layering = infer_volcano_like_layers(&graph, IsogenyGraphNodeId(0));

        assert_eq!(layering.levels, vec![vec![IsogenyGraphNodeId(0)]]);
        assert_eq!(
            role_of(&layering, IsogenyGraphNodeId(0)),
            VolcanoRole::Isolated
        );
    }

    #[test]
    fn simple_single_edge_example_marks_reachable_nodes_as_floor() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let layering = infer_volcano_like_layers(&graph, IsogenyGraphNodeId(0));

        assert_eq!(
            layering.levels,
            vec![vec![IsogenyGraphNodeId(0)], vec![IsogenyGraphNodeId(1)]]
        );
        assert_eq!(
            role_of(&layering, IsogenyGraphNodeId(0)),
            VolcanoRole::Floor
        );
        assert_eq!(
            role_of(&layering, IsogenyGraphNodeId(1)),
            VolcanoRole::Floor
        );
    }

    #[test]
    fn f17_example_exhibits_surface_middle_and_floor_roles() {
        let (graph, root) = graph_with_surface_middle_and_floor_over_f17();
        let layering = infer_volcano_like_layers(&graph, root);
        let roles = layering
            .roles
            .iter()
            .map(|(_, role)| *role)
            .collect::<Vec<_>>();

        assert!(roles.contains(&VolcanoRole::Surface), "roles: {:?}", roles);
        assert!(roles.contains(&VolcanoRole::Middle), "roles: {:?}", roles);
        assert!(roles.contains(&VolcanoRole::Floor), "roles: {:?}", roles);
    }
}
