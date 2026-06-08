use std::fmt;
use std::hash::Hash;

use crate::elliptic_curves::endomorphisms::{
    EndomorphismRingCandidateSet, VolcanoEndomorphismLevelCandidate,
};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId, VolcanoLikeLayering,
    infer_volcano_like_layers,
};
use crate::numerics::PositivePrimeError;
use num_bigint::BigUint;

/// Comparison between arithmetic local-order candidates and the current
/// graph-theoretic volcano heuristic.
///
/// This enum stays intentionally modest. It compares only the number of
/// possible arithmetic local levels with the number of weak-BFS heuristic
/// layers, without pretending that the graph heuristic certifies the true
/// arithmetic volcanic level of the curve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VolcanoHeuristicComparison {
    /// The heuristic produced no usable layers, so no comparison is available.
    HeuristicUnavailable,
    /// The arithmetic candidate levels and the heuristic layering have the
    /// same cardinality.
    CompatibleLevelCount,
    /// The arithmetic candidate levels and the heuristic layering do not have
    /// the same cardinality, so the current comparison remains inconclusive.
    InconclusiveLevelCountMismatch,
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
pub struct EndomorphismVolcanoReport {
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
    pub fn from_graph_heuristic(
        candidate_set: &EndomorphismRingCandidateSet,
        prime: &BigUint,
        graph_heuristic: VolcanoLikeLayering,
    ) -> Result<Self, PositivePrimeError> {
        let local_order_candidates = candidate_set.volcanic_level_candidates_at(prime)?;
        let possible_levels = distinct_levels(&local_order_candidates);
        let comparison_with_graph_heuristic =
            compare_with_heuristic(&local_order_candidates, &graph_heuristic);

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
    pub fn from_graph_and_root<C: GraphCurveModel>(
        candidate_set: &EndomorphismRingCandidateSet,
        prime: &BigUint,
        graph: &IsogenyGraph<C>,
        root: IsogenyGraphNodeId,
    ) -> Result<Self, PositivePrimeError>
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + fmt::Debug,
    {
        let graph_heuristic = infer_volcano_like_layers(graph, root);
        Self::from_graph_heuristic(candidate_set, prime, graph_heuristic)
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the arithmetic local-order candidates `O_f` annotated by `v_ℓ(f)`.
    pub fn local_order_candidates(&self) -> &[VolcanoEndomorphismLevelCandidate] {
        &self.local_order_candidates
    }

    /// Returns the distinct arithmetic levels compatible with the current Frobenius data.
    pub fn possible_levels(&self) -> &[u32] {
        &self.possible_levels
    }

    /// Returns the graph-theoretic weak-BFS volcano heuristic.
    pub fn graph_heuristic(&self) -> &VolcanoLikeLayering {
        &self.graph_heuristic
    }

    /// Returns the modest comparison between arithmetic candidates and graph heuristic.
    pub fn comparison_with_graph_heuristic(&self) -> &VolcanoHeuristicComparison {
        &self.comparison_with_graph_heuristic
    }

    /// Returns the number of heuristic weak-BFS levels.
    pub fn heuristic_level_count(&self) -> usize {
        self.graph_heuristic.levels.len()
    }
}

impl EndomorphismRingCandidateSet {
    /// Builds the bridge report from the current graph-theoretic heuristic.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn volcano_report_with_heuristic(
        &self,
        prime: &BigUint,
        graph_heuristic: VolcanoLikeLayering,
    ) -> Result<EndomorphismVolcanoReport, PositivePrimeError> {
        EndomorphismVolcanoReport::from_graph_heuristic(self, prime, graph_heuristic)
    }

    /// Builds the bridge report by first running the weak-BFS volcano heuristic.
    ///
    /// Complexity: dominated by `num-prime` for the arithmetic candidate annotation
    /// plus the current `Θ(|V| + |E|)` weak-graph traversal performed by
    /// [`infer_volcano_like_layers`].
    pub fn volcano_report_from_graph_and_root<C: GraphCurveModel>(
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

fn distinct_levels(candidates: &[VolcanoEndomorphismLevelCandidate]) -> Vec<u32> {
    let mut levels = candidates
        .iter()
        .map(|candidate| candidate.level())
        .collect::<Vec<_>>();
    levels.sort_unstable();
    levels.dedup();
    levels
}

fn compare_with_heuristic(
    local_order_candidates: &[VolcanoEndomorphismLevelCandidate],
    graph_heuristic: &VolcanoLikeLayering,
) -> VolcanoHeuristicComparison {
    let possible_level_count = distinct_levels(local_order_candidates).len();

    if graph_heuristic.levels.is_empty() {
        return VolcanoHeuristicComparison::HeuristicUnavailable;
    }

    let heuristic_level_count = graph_heuristic.levels.len();
    if heuristic_level_count == possible_level_count {
        VolcanoHeuristicComparison::CompatibleLevelCount
    } else {
        VolcanoHeuristicComparison::InconclusiveLevelCountMismatch
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::QuadraticDiscriminant;
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::{
        EndomorphismVolcanoReport, IsogenyGraphBuilder, IsogenyGraphNodeId,
        VolcanoHeuristicComparison, VolcanoLikeLayering,
    };
    use num_bigint::BigUint;

    type F41 = Fp<41>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn report_can_compare_candidate_levels_with_a_precomputed_heuristic() {
        let candidate_set = QuadraticDiscriminant::new(-16)
            .factorization()
            .expect("-16 should factor canonically")
            .endomorphism_ring_candidates()
            .expect("candidate orders should construct");
        let heuristic = VolcanoLikeLayering {
            levels: vec![vec![IsogenyGraphNodeId(0)], vec![IsogenyGraphNodeId(1)]],
            roles: Vec::new(),
        };

        let report = EndomorphismVolcanoReport::from_graph_heuristic(
            &candidate_set,
            &BigUint::from(2u8),
            heuristic,
        )
        .expect("report should build");

        assert_eq!(report.possible_levels(), &[0, 1]);
        assert_eq!(report.heuristic_level_count(), 2);
        assert_eq!(report.local_order_candidates().len(), 2);
        assert_eq!(
            report.comparison_with_graph_heuristic(),
            &VolcanoHeuristicComparison::CompatibleLevelCount
        );
    }

    #[test]
    fn report_marks_empty_heuristics_as_unavailable() {
        let candidate_set = QuadraticDiscriminant::new(-16)
            .factorization()
            .expect("-16 should factor canonically")
            .endomorphism_ring_candidates()
            .expect("candidate orders should construct");

        let report = EndomorphismVolcanoReport::from_graph_heuristic(
            &candidate_set,
            &BigUint::from(2u8),
            VolcanoLikeLayering {
                levels: Vec::new(),
                roles: Vec::new(),
            },
        )
        .expect("report should build");

        assert_eq!(
            report.comparison_with_graph_heuristic(),
            &VolcanoHeuristicComparison::HeuristicUnavailable
        );
    }

    #[test]
    fn report_can_run_the_existing_graph_heuristic_itself() {
        let candidate_set = QuadraticDiscriminant::new(-16)
            .factorization()
            .expect("-16 should factor canonically")
            .endomorphism_ring_candidates()
            .expect("candidate orders should construct");
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let report = candidate_set
            .volcano_report_from_graph_and_root(&BigUint::from(2u8), &graph, IsogenyGraphNodeId(0))
            .expect("report should build");

        assert_eq!(report.graph_heuristic().levels.len(), 2);
        assert_eq!(
            report.comparison_with_graph_heuristic(),
            &VolcanoHeuristicComparison::CompatibleLevelCount
        );
    }
}
