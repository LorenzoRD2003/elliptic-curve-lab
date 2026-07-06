use std::collections::BTreeSet;
use std::hash::Hash;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId, VolcanoRole,
    endomorphisms::IsogenyEdgeEndomorphismTentativeRelation,
};

/// Narrow graph-side `ℓ`-volcano evidence extracted while building an
/// endomorphism report.
///
/// The arithmetic candidate sets remain the source of truth. This helper only
/// records endpoint-role constraints and surface-anchored edge directions that
/// are useful for educational refinement; it does not certify the true
/// volcanic level of any node.
pub(super) struct ObservedGraphVolcanoEvidence {
    levels: Vec<Option<u32>>,
    roles: Vec<VolcanoRole>,
    surface_anchored: bool,
}

impl ObservedGraphVolcanoEvidence {
    pub(super) fn from_graph<C: GraphCurveModel>(graph: &IsogenyGraph<C>) -> Self
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + core::fmt::Debug,
    {
        let root = IsogenyGraphNodeId(0);
        let layering = graph.infer_volcano_like_layers(root);
        let mut levels = vec![None; graph.node_count()];
        for (level, node_ids) in layering.levels().iter().enumerate() {
            for node_id in node_ids {
                levels[node_id.0] = Some(level as u32);
            }
        }
        let mut roles = vec![VolcanoRole::Unknown; graph.node_count()];
        for (node_id, role) in layering.roles() {
            roles[node_id.0] = *role;
        }
        let surface_anchored = layering.role_of(root) == Some(VolcanoRole::Surface);

        Self {
            levels,
            roles,
            surface_anchored,
        }
    }

    pub(super) fn allowed_levels_for(
        &self,
        node_id: IsogenyGraphNodeId,
        possible_levels: &[u32],
    ) -> Option<BTreeSet<u32>> {
        match self.roles.get(node_id.0).copied()? {
            VolcanoRole::Surface if possible_levels.contains(&0) => Some(BTreeSet::from([0])),
            VolcanoRole::Floor => possible_levels
                .iter()
                .copied()
                .max()
                .map(|max_level| BTreeSet::from([max_level])),
            VolcanoRole::Middle | VolcanoRole::Isolated | VolcanoRole::Unknown => None,
            VolcanoRole::Surface => None,
        }
    }

    pub(super) fn relation_for(
        &self,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
    ) -> Option<IsogenyEdgeEndomorphismTentativeRelation> {
        if !self.surface_anchored {
            return None;
        }

        let source_level = self.levels.get(source.0).copied().flatten()?;
        let target_level = self.levels.get(target.0).copied().flatten()?;
        if source_level == target_level {
            Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal)
        } else if source_level == target_level + 1 {
            Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending)
        } else if target_level == source_level + 1 {
            Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending)
        } else {
            None
        }
    }
}
