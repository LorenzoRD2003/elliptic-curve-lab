use num_bigint::BigUint;
use std::fmt;
use std::hash::Hash;

use crate::elliptic_curves::endomorphisms::candidate_sets::{
    EndomorphismRingCandidateSet, VolcanoEndomorphismLevelCandidate,
};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId, VolcanoLikeLayering,
};
use crate::numerics::PositivePrimeError;

/// Comparison between arithmetic local-order candidates and the current
/// graph-theoretic volcano heuristic.
///
/// This enum stays intentionally modest. It compares only the number of
/// possible arithmetic local levels with the number of weak-BFS heuristic
/// layers, without pretending that the graph heuristic certifies the true
/// arithmetic volcanic level of the curve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum VolcanoHeuristicComparison {
    /// The heuristic produced no usable layers, so no comparison is available.
    HeuristicUnavailable,
    /// The arithmetic candidate levels and the heuristic layering have the
    /// same cardinality.
    CompatibleLevelCount,
    /// The arithmetic candidate levels and the heuristic layering do not have
    /// the same cardinality, so the current comparison remains inconclusive.
    InconclusiveLevelCountMismatch,
}

impl VolcanoHeuristicComparison {
    fn from_evidence(
        local_order_candidates: &[VolcanoEndomorphismLevelCandidate],
        graph_heuristic: &VolcanoLikeLayering,
    ) -> Self {
        let possible_level_count =
            VolcanoEndomorphismLevelCandidate::distinct_levels_from(local_order_candidates).len();

        if graph_heuristic.is_empty() {
            return Self::HeuristicUnavailable;
        }

        let heuristic_level_count = graph_heuristic.level_count();
        if heuristic_level_count == possible_level_count {
            Self::CompatibleLevelCount
        } else {
            Self::InconclusiveLevelCountMismatch
        }
    }
}

/// Bridge report between arithmetic endomorphism-order candidates and the
/// current graph-theoretic volcano heuristic.
///
/// This report does not replace [`infer_volcano_like_layers`]. Instead, it
/// packages:
///
/// - the `ℓ`-local order candidates derived from Frobenius-side arithmetic
/// - the distinct possible arithmetic levels `v_ℓ(f)`
/// - one weak-BFS graph heuristic layering
/// - a deliberately modest comparison between those two layers of evidence
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct EndomorphismVolcanoReport {
    prime: BigUint,
    local_order_candidates: Vec<VolcanoEndomorphismLevelCandidate>,
    possible_levels: Vec<u32>,
    graph_heuristic: VolcanoLikeLayering,
    comparison_with_graph_heuristic: VolcanoHeuristicComparison,
}

impl EndomorphismVolcanoReport {
    /// Builds the report from an already computed graph-theoretic heuristic.
    ///
    /// Complexity: dominated by `num-prime`.
    pub(crate) fn from_graph_heuristic(
        candidate_set: &EndomorphismRingCandidateSet,
        prime: &BigUint,
        graph_heuristic: VolcanoLikeLayering,
    ) -> Result<Self, PositivePrimeError> {
        let local_order_candidates = candidate_set.volcanic_level_candidates_at(prime)?;
        let possible_levels =
            VolcanoEndomorphismLevelCandidate::distinct_levels_from(&local_order_candidates);
        let comparison_with_graph_heuristic =
            VolcanoHeuristicComparison::from_evidence(&local_order_candidates, &graph_heuristic);

        Ok(Self {
            prime: prime.clone(),
            local_order_candidates,
            possible_levels,
            graph_heuristic,
            comparison_with_graph_heuristic,
        })
    }

    /// Builds the report by first running the current weak-BFS volcano heuristic.
    ///
    /// Complexity: dominated by `num-prime` for the arithmetic candidate annotation
    /// plus the current `Θ(|V| + |E|)` weak-graph traversal performed by
    /// [`infer_volcano_like_layers`].
    pub(crate) fn from_graph_and_root<C: GraphCurveModel>(
        candidate_set: &EndomorphismRingCandidateSet,
        prime: &BigUint,
        graph: &IsogenyGraph<C>,
        root: IsogenyGraphNodeId,
    ) -> Result<Self, PositivePrimeError>
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + fmt::Debug,
    {
        let graph_heuristic = graph.infer_volcano_like_layers(root);
        Self::from_graph_heuristic(candidate_set, prime, graph_heuristic)
    }

    /// Returns how many arithmetic local-order candidates `O_f` were recorded.
    pub(crate) fn local_order_candidate_count(&self) -> usize {
        self.local_order_candidates.len()
    }

    /// Returns the distinct arithmetic levels compatible with the current Frobenius data.
    pub(crate) fn possible_levels(&self) -> &[u32] {
        &self.possible_levels
    }

    /// Returns the graph-theoretic weak-BFS volcano heuristic.
    pub(crate) fn graph_heuristic(&self) -> &VolcanoLikeLayering {
        &self.graph_heuristic
    }

    /// Returns the modest comparison between arithmetic candidates and graph heuristic.
    pub(crate) fn comparison_with_graph_heuristic(&self) -> &VolcanoHeuristicComparison {
        &self.comparison_with_graph_heuristic
    }

    /// Returns the number of heuristic weak-BFS levels.
    pub(crate) fn heuristic_level_count(&self) -> usize {
        self.graph_heuristic.level_count()
    }
}

impl EndomorphismRingCandidateSet {
    /// Builds the bridge report by first running the weak-BFS volcano heuristic.
    ///
    /// Complexity: dominated by `num-prime` for the arithmetic candidate annotation
    /// plus the current `Θ(|V| + |E|)` weak-graph traversal performed by
    /// [`infer_volcano_like_layers`].
    pub(crate) fn volcano_report_from_graph_and_root<C: GraphCurveModel>(
        &self,
        prime: &BigUint,
        graph: &IsogenyGraph<C>,
        root: IsogenyGraphNodeId,
    ) -> Result<EndomorphismVolcanoReport, PositivePrimeError>
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + fmt::Debug,
    {
        EndomorphismVolcanoReport::from_graph_and_root(self, prime, graph, root)
    }
}

#[cfg(test)]
mod tests;
