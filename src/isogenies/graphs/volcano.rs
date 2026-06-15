use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId};

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
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VolcanoLikeLayering {
    levels: Vec<Vec<IsogenyGraphNodeId>>,
    roles: Vec<(IsogenyGraphNodeId, VolcanoRole)>,
}

impl VolcanoLikeLayering {
    pub(crate) fn new(
        levels: Vec<Vec<IsogenyGraphNodeId>>,
        roles: Vec<(IsogenyGraphNodeId, VolcanoRole)>,
    ) -> Self {
        Self { levels, roles }
    }

    /// Returns the weak-BFS levels in distance order from the chosen root.
    pub fn levels(&self) -> &[Vec<IsogenyGraphNodeId>] {
        &self.levels
    }

    /// Returns the stored `(node, role)` annotations.
    pub fn roles(&self) -> &[(IsogenyGraphNodeId, VolcanoRole)] {
        &self.roles
    }

    /// Returns the nodes recorded at one weak-BFS level, if that level exists.
    pub fn nodes_at_level(&self, level: usize) -> Option<&[IsogenyGraphNodeId]> {
        self.levels.get(level).map(Vec::as_slice)
    }

    /// Returns the recorded role of one node when present in the layering data.
    pub fn role_of(&self, node: IsogenyGraphNodeId) -> Option<VolcanoRole> {
        self.roles
            .iter()
            .find_map(|(candidate, role)| (*candidate == node).then_some(*role))
    }

    /// Counts how many nodes currently carry the requested role.
    pub fn count_role(&self, role: VolcanoRole) -> usize {
        self.roles
            .iter()
            .filter(|(_, candidate)| *candidate == role)
            .count()
    }

    /// Returns how many weak-BFS levels were recorded.
    pub fn level_count(&self) -> usize {
        self.levels.len()
    }

    /// Returns whether the heuristic produced no reachable layering.
    pub fn is_empty(&self) -> bool {
        self.levels.is_empty()
    }
}

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
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
    pub fn infer_volcano_like_layers(&self, root: IsogenyGraphNodeId) -> VolcanoLikeLayering {
        let Some(_) = self.node(root) else {
            return VolcanoLikeLayering::new(
                Vec::new(),
                self.nodes()
                    .iter()
                    .map(|node| (node.id(), VolcanoRole::Unknown))
                    .collect(),
            );
        };

        let distances = self.weak_bfs_distances(root);
        let levels = levels_from_distances(&distances);
        let surface_level = self
            .nodes()
            .iter()
            .filter_map(|node| {
                let id = node.id();
                distances[id.0].filter(|_| self.is_high_degree_node(id))
            })
            .min();

        let roles = self
            .nodes()
            .iter()
            .map(|node| {
                let id = node.id();
                let role = match distances[id.0] {
                    None => VolcanoRole::Unknown,
                    Some(level) => self.classify_role(id, level, surface_level),
                };
                (id, role)
            })
            .collect();

        VolcanoLikeLayering::new(levels, roles)
    }

    fn classify_role(
        &self,
        node: IsogenyGraphNodeId,
        level: usize,
        surface_level: Option<usize>,
    ) -> VolcanoRole {
        let total = self.weak_degree(node);
        if total == 0 {
            VolcanoRole::Isolated
        } else if total == 1 && !self.has_self_loop(node) {
            VolcanoRole::Floor
        } else if surface_level == Some(level) && self.is_high_degree_node(node) {
            VolcanoRole::Surface
        } else {
            VolcanoRole::Middle
        }
    }
    fn is_high_degree_node(&self, node: IsogenyGraphNodeId) -> bool {
        self.weak_degree(node) > 1 || self.has_self_loop(node)
    }
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

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Fp, traits::Field};
    use crate::isogenies::graphs::{
        IsogenyGraph, IsogenyGraphBuilder, IsogenyGraphNodeId, VolcanoRole,
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
                    let roles = graph
                        .infer_volcano_like_layers(root)
                        .roles()
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
        layering: &crate::isogenies::graphs::volcano::VolcanoLikeLayering,
        node: IsogenyGraphNodeId,
    ) -> VolcanoRole {
        layering
            .role_of(node)
            .expect("test node should have a role")
    }

    #[test]
    fn invalid_root_marks_all_nodes_unknown() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let layering = graph.infer_volcano_like_layers(IsogenyGraphNodeId(99));

        assert!(layering.is_empty());
        assert!(
            layering
                .roles()
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

        let layering = graph.infer_volcano_like_layers(IsogenyGraphNodeId(0));

        assert_eq!(layering.levels(), &[vec![IsogenyGraphNodeId(0)]]);
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

        let layering = graph.infer_volcano_like_layers(IsogenyGraphNodeId(0));

        assert_eq!(
            layering.levels(),
            &[vec![IsogenyGraphNodeId(0)], vec![IsogenyGraphNodeId(1)]]
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
        let layering = graph.infer_volcano_like_layers(root);
        let roles = layering
            .roles()
            .iter()
            .map(|(_, role)| *role)
            .collect::<Vec<_>>();

        assert!(roles.contains(&VolcanoRole::Surface), "roles: {:?}", roles);
        assert!(roles.contains(&VolcanoRole::Middle), "roles: {:?}", roles);
        assert!(roles.contains(&VolcanoRole::Floor), "roles: {:?}", roles);
    }
}
