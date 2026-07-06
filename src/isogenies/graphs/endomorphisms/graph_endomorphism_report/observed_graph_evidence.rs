use std::collections::BTreeSet;
use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyEdgeEndomorphismTentativeRelation, floor_search::VolcanoAltimeterEvidence,
        graph_endomorphism_report::observed_volcano_evidence::ObservedGraphVolcanoEvidence,
    },
};

/// Combined optional graph-side evidence used while building endomorphism reports.
///
/// Altimeter evidence has priority when both endpoint floor distances `δ(v)`
/// are certified. The older weak-BFS volcano evidence remains as a fallback for
/// educational cases where the altimeter cannot certify both endpoints.
pub(super) struct ObservedEndomorphismGraphEvidence {
    weak_volcano: ObservedGraphVolcanoEvidence,
    altimeter: VolcanoAltimeterEvidence,
}

impl ObservedEndomorphismGraphEvidence {
    pub(super) fn from_graph<C: GraphCurveModel>(graph: &IsogenyGraph<C>, prime: &BigUint) -> Self
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + fmt::Debug,
    {
        Self {
            weak_volcano: ObservedGraphVolcanoEvidence::from_graph(graph),
            altimeter: VolcanoAltimeterEvidence::from_graph(graph, prime),
        }
    }

    pub(super) fn allowed_levels_for(
        &self,
        node_id: IsogenyGraphNodeId,
        possible_levels: &[u32],
    ) -> Option<BTreeSet<u32>> {
        self.weak_volcano
            .allowed_levels_for(node_id, possible_levels)
    }

    pub(super) fn edge_relation_for(
        &self,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
        arithmetic_relation: &IsogenyEdgeEndomorphismTentativeRelation,
    ) -> Option<IsogenyEdgeEndomorphismTentativeRelation> {
        if arithmetic_relation == &IsogenyEdgeEndomorphismTentativeRelation::Unsupported {
            return None;
        }

        self.altimeter
            .relation_for(source, target)
            .or_else(|| self.weak_volcano.relation_for(source, target))
    }
}
