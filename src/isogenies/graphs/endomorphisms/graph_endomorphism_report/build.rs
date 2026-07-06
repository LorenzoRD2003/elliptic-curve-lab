use std::hash::Hash;

use num_bigint::BigUint;

use crate::elliptic_curves::{
    endomorphisms::candidate_sets::VolcanoEndomorphismLevelCandidate,
    frobenius::FrobeniusTraceCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphError,
    endomorphisms::graph_endomorphism_report::{
        IsogenyGraphEndomorphismEdgeReport, IsogenyGraphEndomorphismNodeReport,
        IsogenyGraphEndomorphismReport, observed_graph_evidence::ObservedEndomorphismGraphEvidence,
    },
};

impl IsogenyGraphEndomorphismReport {
    /// Builds the report from one graph and one chosen prime `ℓ`.
    ///
    /// Complexity:
    /// - one exhaustive Frobenius-trace computation per node
    /// - arithmetic dominated by `num-prime` for each node
    /// - one tentative edge comparison per stored edge
    pub(crate) fn from_graph<C: GraphCurveModel + FrobeniusTraceCurveModel>(
        graph: &IsogenyGraph<C>,
        prime: &BigUint,
    ) -> Result<Self, IsogenyGraphError>
    where
        C::BaseField:
            EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
        C::Point: Clone + Eq + Hash + PartialEq,
        C::IsomorphismWitness: Clone + core::fmt::Debug,
    {
        let node_candidate_sets = graph.graph_endomorphism_candidates()?;

        let graph_evidence = ObservedEndomorphismGraphEvidence::from_graph(graph, prime);

        let nodes = node_candidate_sets
            .iter()
            .map(|(node_id, candidate_set)| {
                let local_levels = candidate_set
                    .volcanic_level_candidates_at(prime)
                    .map_err(|_| IsogenyGraphError::InvalidDegree)?;
                let possible_levels =
                    VolcanoEndomorphismLevelCandidate::distinct_levels_from(&local_levels);
                let observed_allowed_levels =
                    graph_evidence.allowed_levels_for(*node_id, &possible_levels);
                Ok(IsogenyGraphEndomorphismNodeReport::new(
                    *node_id,
                    candidate_set.clone(),
                    local_levels,
                    observed_allowed_levels,
                ))
            })
            .collect::<Result<Vec<_>, IsogenyGraphError>>()?;

        let edges = graph
            .edges()
            .iter()
            .map(|edge| {
                let source = edge.source();
                let target = edge.target();
                let source_candidates = nodes[source.0].candidate_set();
                let target_candidates = nodes[target.0].candidate_set();
                let relation = source_candidates
                    .tentative_edge_endomorphism_report(prime, target_candidates)
                    .map_err(|_| IsogenyGraphError::InvalidDegree)?;
                let observed_relation =
                    graph_evidence.edge_relation_for(source, target, relation.relation());

                Ok(IsogenyGraphEndomorphismEdgeReport::new(
                    edge.id(),
                    source,
                    target,
                    relation,
                    observed_relation,
                ))
            })
            .collect::<Result<Vec<_>, IsogenyGraphError>>()?;

        Ok(Self::new(prime.clone(), nodes, edges))
    }
}

impl<C: GraphCurveModel + FrobeniusTraceCurveModel> IsogenyGraph<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + Hash + PartialEq,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    /// Builds the educational endomorphism-side report for this `ℓ`-isogeny graph.
    ///
    /// The report annotates every stored node with the Frobenius-compatible
    /// imaginary quadratic orders currently possible for that representative,
    /// then compares source and target `ℓ`-local conductor levels along each
    /// stored edge.
    ///
    /// This is intentionally a tentative report: it does not certify the exact
    /// endomorphism ring of any curve, and it does not prove definitive
    /// horizontal/ascending/descending edge types.
    ///
    /// Complexity:
    /// - one exhaustive Frobenius-trace computation per node
    /// - arithmetic dominated by `num-prime` for each node
    /// - one tentative edge comparison per stored edge
    pub fn endomorphism_report_at(
        &self,
        prime: &BigUint,
    ) -> Result<IsogenyGraphEndomorphismReport, IsogenyGraphError> {
        IsogenyGraphEndomorphismReport::from_graph(self, prime)
    }
}
