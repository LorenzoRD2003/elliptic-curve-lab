use num_bigint::BigUint;

use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{
        ShortestFloorPathReport,
        volcano_structure::{
            UncertifiedVolcanoNodeReport, VolcanoStructureLevelReport, VolcanoStructureNodeReport,
            VolcanoStructureRole, VolcanoStructureUncertifiedReason,
        },
    },
};

/// Structural report for a stored graph viewed as an ordinary `ℓ`-volcano.
///
/// The report is built from certified shortest paths to the floor. For each
/// certified node it records `δ(v) = dist(v, V_d)`. The largest certified
/// distance is used as the certified depth `d̂`, and the node's certified level
/// is then `d̂ - δ(v)`.
///
/// This is a structural report for the graph evidence that is actually stored,
/// not a proof that the ambient isogeny component has been completely
/// enumerated. Nodes whose floor distance cannot be certified are kept in
/// [`Self::uncertified_nodes`] rather than forced into a level.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VolcanoStructureReport {
    prime: BigUint,
    certified_depth: Option<usize>,
    levels: Vec<VolcanoStructureLevelReport>,
    certified_nodes: Vec<VolcanoStructureNodeReport>,
    uncertified_nodes: Vec<UncertifiedVolcanoNodeReport>,
}

impl VolcanoStructureReport {
    pub(crate) fn from_floor_paths(
        prime: BigUint,
        floor_paths: Vec<ShortestFloorPathReport>,
        uncertified_nodes: Vec<UncertifiedVolcanoNodeReport>,
    ) -> Self {
        let certified_depth = floor_paths
            .iter()
            .map(ShortestFloorPathReport::distance_to_floor)
            .max();

        let mut certified_nodes = floor_paths
            .into_iter()
            .map(|floor_path| {
                VolcanoStructureNodeReport::from_floor_path(
                    certified_depth.expect("at least one path exists while mapping paths"),
                    floor_path,
                )
            })
            .collect::<Vec<_>>();
        certified_nodes.sort_by_key(|node| node.node_id().0);

        let levels = match certified_depth {
            Some(depth) => Self::levels_from_certified_nodes(depth, &certified_nodes),
            None => Vec::new(),
        };

        Self {
            prime,
            certified_depth,
            levels,
            certified_nodes,
            uncertified_nodes,
        }
    }

    /// Returns the chosen local prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the maximum depth `d̂` certified by stored floor-distance evidence.
    ///
    /// This is the largest certified `δ(v)` among stored nodes. It coincides
    /// with the volcano depth `d` when the stored graph evidence includes a
    /// certified surface vertex in a complete ordinary component.
    pub fn certified_depth(&self) -> Option<usize> {
        self.certified_depth
    }

    /// Returns the certified levels in ascending order from surface to floor.
    pub fn levels(&self) -> &[VolcanoStructureLevelReport] {
        &self.levels
    }

    /// Returns the certified surface `V₀`, also called the crater.
    pub fn surface(&self) -> Option<&VolcanoStructureLevelReport> {
        self.level(0)
    }

    /// Returns the certified crater, equal to the surface `V₀`.
    pub fn crater(&self) -> Option<&VolcanoStructureLevelReport> {
        self.surface()
    }

    /// Returns the deepest certified level `V_d`.
    pub fn floor(&self) -> Option<&VolcanoStructureLevelReport> {
        self.certified_depth.and_then(|depth| self.level(depth))
    }

    /// Returns one certified level by index.
    pub fn level(&self, level: usize) -> Option<&VolcanoStructureLevelReport> {
        self.levels
            .get(level)
            .filter(|report| report.level() == level)
    }

    /// Returns the certified node reports in dense node-id order.
    pub fn certified_nodes(&self) -> &[VolcanoStructureNodeReport] {
        &self.certified_nodes
    }

    /// Returns the nodes whose floor distance could not be certified.
    pub fn uncertified_nodes(&self) -> &[UncertifiedVolcanoNodeReport] {
        &self.uncertified_nodes
    }

    /// Returns whether every stored node has certified volcano-structure data.
    pub fn is_fully_certified(&self) -> bool {
        self.uncertified_nodes.is_empty()
    }

    /// Returns one certified node report when available.
    pub fn certified_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Option<&VolcanoStructureNodeReport> {
        self.certified_nodes
            .iter()
            .find(|node| node.node_id() == node_id)
    }

    /// Returns one uncertified node report when available.
    pub fn uncertified_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Option<&UncertifiedVolcanoNodeReport> {
        self.uncertified_nodes
            .iter()
            .find(|node| node.node_id() == node_id)
    }

    pub(crate) fn partial_uncertified_node_ids(&self) -> Vec<IsogenyGraphNodeId> {
        self.uncertified_nodes
            .iter()
            .filter_map(|node| {
                (node.reason() == &VolcanoStructureUncertifiedReason::PartialGraph)
                    .then_some(node.node_id())
            })
            .collect()
    }

    fn levels_from_certified_nodes(
        certified_depth: usize,
        certified_nodes: &[VolcanoStructureNodeReport],
    ) -> Vec<VolcanoStructureLevelReport> {
        let mut level_nodes = vec![Vec::new(); certified_depth + 1];
        for node in certified_nodes {
            level_nodes[node.level()].push(node.node_id());
        }

        level_nodes
            .into_iter()
            .enumerate()
            .map(|(level, nodes)| {
                VolcanoStructureLevelReport::new(
                    level,
                    VolcanoStructureRole::for_level(level, certified_depth),
                    nodes,
                )
            })
            .collect()
    }
}
