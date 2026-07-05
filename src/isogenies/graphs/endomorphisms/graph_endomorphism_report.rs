use num_bigint::BigUint;
use std::hash::Hash;

use crate::elliptic_curves::{
    endomorphisms::candidate_sets::{
        EndomorphismRingCandidateSet, VolcanoEndomorphismLevelCandidate,
    },
    frobenius::FrobeniusTraceCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphEdgeId, IsogenyGraphError, IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyEdgeEndomorphismReport,
        refinement::{
            CandidateRefinementError, CandidateRefinementStrategy, EndomorphismCandidateRefinement,
            IncidentEdgeRefinementConstraint,
        },
    },
};

/// Endomorphism-side report for one stored graph node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismNodeReport {
    node_id: IsogenyGraphNodeId,
    candidate_set: EndomorphismRingCandidateSet,
    local_levels: Vec<VolcanoEndomorphismLevelCandidate>,
}

/// Tentative endomorphism-side report for one stored graph edge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismEdgeReport {
    edge_id: IsogenyGraphEdgeId,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    relation: IsogenyEdgeEndomorphismReport,
}

/// Endomorphism-side report for an entire educational isogeny graph at one chosen prime `ℓ`.
///
/// This aggregate report is still conservative. It packages:
///
/// - automatic Frobenius-compatible candidate-order data for each node
/// - the corresponding `ℓ`-local candidate levels at each node
/// - tentative edge relations derived from those node-wise candidate sets
///
/// It does **not** certify exact endomorphism rings or definitive edge types.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismReport {
    prime: BigUint,
    nodes: Vec<IsogenyGraphEndomorphismNodeReport>,
    edges: Vec<IsogenyGraphEndomorphismEdgeReport>,
}

impl IsogenyGraphEndomorphismNodeReport {
    pub(crate) fn new(
        node_id: IsogenyGraphNodeId,
        candidate_set: EndomorphismRingCandidateSet,
        local_levels: Vec<VolcanoEndomorphismLevelCandidate>,
    ) -> Self {
        Self {
            node_id,
            candidate_set,
            local_levels,
        }
    }

    /// Returns the node identifier.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the Frobenius-compatible candidate orders for this node.
    pub fn candidate_set(&self) -> &EndomorphismRingCandidateSet {
        &self.candidate_set
    }

    /// Returns how many arithmetic `ℓ`-local level candidates were recorded
    /// for this node.
    pub fn local_level_candidate_count(&self) -> usize {
        self.local_levels.len()
    }

    /// Returns the distinct possible local levels for this node.
    pub fn possible_levels(&self) -> Vec<u32> {
        VolcanoEndomorphismLevelCandidate::distinct_levels_from(&self.local_levels)
    }
}

impl IsogenyGraphEndomorphismEdgeReport {
    /// Returns the edge identifier.
    pub fn edge_id(&self) -> IsogenyGraphEdgeId {
        self.edge_id
    }

    /// Returns the source node identifier.
    pub fn source(&self) -> IsogenyGraphNodeId {
        self.source
    }

    /// Returns the target node identifier.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the tentative endomorphism-side edge relation report.
    pub fn relation(&self) -> &IsogenyEdgeEndomorphismReport {
        &self.relation
    }
}

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

        let nodes = node_candidate_sets
            .iter()
            .map(|(node_id, candidate_set)| {
                let local_levels = candidate_set
                    .volcanic_level_candidates_at(prime)
                    .map_err(|_| IsogenyGraphError::InvalidDegree)?;
                Ok(IsogenyGraphEndomorphismNodeReport::new(
                    *node_id,
                    candidate_set.clone(),
                    local_levels,
                ))
            })
            .collect::<Result<Vec<_>, IsogenyGraphError>>()?;

        let edges = graph
            .edges()
            .iter()
            .map(|edge| {
                let source = edge.source();
                let target = edge.target();
                let source_candidates = &nodes[source.0].candidate_set;
                let target_candidates = &nodes[target.0].candidate_set;
                let relation = source_candidates
                    .tentative_edge_endomorphism_report(prime, target_candidates)
                    .map_err(|_| IsogenyGraphError::InvalidDegree)?;

                Ok(IsogenyGraphEndomorphismEdgeReport {
                    edge_id: edge.id(),
                    source,
                    target,
                    relation,
                })
            })
            .collect::<Result<Vec<_>, IsogenyGraphError>>()?;

        Ok(Self {
            prime: prime.clone(),
            nodes,
            edges,
        })
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the node reports in dense node-id order.
    pub fn nodes(&self) -> &[IsogenyGraphEndomorphismNodeReport] {
        &self.nodes
    }

    /// Returns the edge reports in stored edge order.
    pub fn edges(&self) -> &[IsogenyGraphEndomorphismEdgeReport] {
        &self.edges
    }

    /// Returns the node report for the requested id when present.
    pub fn node_report(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Option<&IsogenyGraphEndomorphismNodeReport> {
        self.nodes
            .get(node_id.0)
            .filter(|report| report.node_id == node_id)
    }

    /// Refines the candidate endomorphism orders for one node using evidence
    /// already recorded in this graph report.
    ///
    /// [`CandidateRefinementStrategy::NodeLocalLevelsOnly`] uses only the
    /// node-local conductor-level constraint `v_ℓ(f) ∈ L`.
    /// [`CandidateRefinementStrategy::Conservative`] and
    /// [`CandidateRefinementStrategy::IncidentUnambiguousEdges`] also use
    /// one-hop incident edge constraints, but only when the stored edge
    /// relation is unequivocal. `Ambiguous` and `Unsupported` edge evidence is
    /// ignored by these conservative strategies.
    ///
    /// The result is not a certification of `End(E)`: even a unique survivor is
    /// only the unique candidate compatible with the evidence used by the
    /// selected strategy.
    pub fn refine_candidates_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
        strategy: CandidateRefinementStrategy,
    ) -> Result<EndomorphismCandidateRefinement, CandidateRefinementError> {
        let node_report = self
            .node_report(node_id)
            .ok_or(CandidateRefinementError::NodeNotFound { node_id })?;

        match strategy {
            CandidateRefinementStrategy::NodeLocalLevelsOnly => {
                EndomorphismCandidateRefinement::from_node_local_levels(node_report, self.prime())
            }
            CandidateRefinementStrategy::Conservative
            | CandidateRefinementStrategy::IncidentUnambiguousEdges => {
                EndomorphismCandidateRefinement::from_incident_unambiguous_edges(
                    node_report,
                    self.prime(),
                    self.incident_edge_constraints_for_node(node_id)?,
                )
            }
        }
    }

    fn incident_edge_constraints_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Result<Vec<IncidentEdgeRefinementConstraint>, CandidateRefinementError> {
        let mut constraints = Vec::new();

        for edge in &self.edges {
            let adjacent_node = if edge.source() == node_id {
                edge.target()
            } else if edge.target() == node_id {
                edge.source()
            } else {
                continue;
            };

            let adjacent_levels = self
                .node_report(adjacent_node)
                .ok_or(CandidateRefinementError::NodeNotFound {
                    node_id: adjacent_node,
                })?
                .possible_levels()
                .into_iter()
                .collect();

            if let Some(constraint) =
                IncidentEdgeRefinementConstraint::from_edge_report(node_id, edge, adjacent_levels)
            {
                constraints.push(constraint);
            }
        }

        Ok(constraints)
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

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::isogenies::graphs::{
        IsogenyGraphBuilder, IsogenyGraphNodeId,
        endomorphisms::{
            IsogenyEdgeEndomorphismTentativeRelation, IsogenyGraphEndomorphismReport,
            refinement::{
                CandidateRefinementError, CandidateRefinementStrategy, ConstraintSource,
                EndomorphismCandidateRefinement, LocalEndomorphismConstraint, RefinementConfidence,
            },
        },
    };

    type F41 = crate::fields::Fp41;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn graph_report(
        depth: usize,
    ) -> (
        crate::isogenies::graphs::IsogenyGraph<Curve41>,
        IsogenyGraphEndomorphismReport,
    ) {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(depth)
            .build()
            .expect("graph should build from the concrete curve");
        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");
        (graph, report)
    }

    fn refinement_at(
        depth: usize,
        strategy: CandidateRefinementStrategy,
    ) -> EndomorphismCandidateRefinement {
        let (_, report) = graph_report(depth);
        report
            .refine_candidates_for_node(IsogenyGraphNodeId(0), strategy)
            .expect("candidate refinement should build")
    }

    fn assert_only_node_or_unambiguous_edge_constraints(
        refinement: &EndomorphismCandidateRefinement,
    ) {
        assert!(refinement.constraints().iter().all(|constraint| {
            matches!(constraint, LocalEndomorphismConstraint::NodeLevel { .. })
                || matches!(
                    constraint,
                    LocalEndomorphismConstraint::EdgeRelation {
                        relation: IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
                            | IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending
                            | IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending,
                        ..
                    }
                )
        }));
    }

    fn assert_refinement_candidates_stay_inside_initial_set(
        refinement: &EndomorphismCandidateRefinement,
    ) {
        let initial_orders = refinement.initial_candidates().candidate_orders();
        assert!(!initial_orders.is_empty());
        assert!(
            refinement
                .surviving_candidates()
                .iter()
                .all(|candidate| initial_orders.contains(candidate))
        );
        assert!(
            refinement
                .eliminated_candidates()
                .iter()
                .all(|elimination| initial_orders.contains(elimination.candidate()))
        );
    }

    #[test]
    fn graph_report_collects_node_and_edge_endomorphism_data() {
        let (graph, report) = graph_report(1);

        assert_eq!(report.prime(), &BigUint::from(2u8));
        assert_eq!(report.nodes().len(), graph.node_count());
        assert_eq!(report.edges().len(), graph.edge_count());
        assert!(
            report
                .nodes()
                .iter()
                .all(|node| !node.candidate_set().is_empty())
        );
        assert!(report.nodes().iter().all(
            |node| node.local_level_candidate_count() > 0 && !node.possible_levels().is_empty()
        ));
    }

    #[test]
    fn graph_report_exposes_dense_node_lookup() {
        let (_, report) = graph_report(0);

        let node_report = report
            .node_report(IsogenyGraphNodeId(0))
            .expect("root node report should exist");

        assert_eq!(node_report.node_id(), IsogenyGraphNodeId(0));
        assert_eq!(
            node_report.candidate_set(),
            report.nodes()[0].candidate_set()
        );
        assert!(report.node_report(IsogenyGraphNodeId(99)).is_none());
    }

    #[test]
    fn graph_report_edge_endomorphism_relations_are_tentative() {
        let (_, report) = graph_report(1);

        assert!(report.edges().iter().all(|edge| {
            matches!(
                edge.relation().relation(),
                IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
                    | IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending
                    | IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending
                    | IsogenyEdgeEndomorphismTentativeRelation::Ambiguous
                    | IsogenyEdgeEndomorphismTentativeRelation::Unsupported
            )
        }));
    }

    #[test]
    fn node_local_refinement_keeps_candidates_supported_by_node_levels() {
        let (_, report) = graph_report(0);
        let node_report = report
            .node_report(IsogenyGraphNodeId(0))
            .expect("root node report should exist");

        let refinement = refinement_at(0, CandidateRefinementStrategy::NodeLocalLevelsOnly);

        assert_eq!(refinement.node_id(), IsogenyGraphNodeId(0));
        assert_eq!(refinement.initial_candidates(), node_report.candidate_set());
        assert_eq!(
            refinement.surviving_candidates(),
            node_report.candidate_set().candidate_orders()
        );
        assert!(refinement.eliminated_candidates().is_empty());
        assert_eq!(
            refinement.confidence(),
            RefinementConfidence::ConservativeLocalEvidence
        );
        assert_eq!(refinement.constraints().len(), 1);

        let LocalEndomorphismConstraint::NodeLevel {
            ell,
            allowed_levels,
            provenance,
        } = &refinement.constraints()[0]
        else {
            panic!("node-local refinement should record one node-level constraint");
        };
        let expected_levels = node_report.possible_levels().into_iter().collect();
        assert_eq!(ell, &BigUint::from(2u8));
        assert_eq!(allowed_levels, &expected_levels);
        assert_eq!(
            provenance,
            &ConstraintSource::NodeReport {
                node_id: IsogenyGraphNodeId(0)
            }
        );
    }

    #[test]
    fn conservative_refinement_matches_node_local_evidence_without_incident_edges() {
        let conservative = refinement_at(0, CandidateRefinementStrategy::Conservative);
        let node_local = refinement_at(0, CandidateRefinementStrategy::NodeLocalLevelsOnly);

        assert_eq!(conservative, node_local);
    }

    #[test]
    fn incident_edge_refinement_uses_only_unambiguous_edge_constraints() {
        let refinement = refinement_at(1, CandidateRefinementStrategy::IncidentUnambiguousEdges);

        assert_only_node_or_unambiguous_edge_constraints(&refinement);
    }

    #[test]
    fn conservative_refinement_runs_end_to_end_from_curve_graph() {
        let refinement = refinement_at(1, CandidateRefinementStrategy::Conservative);

        assert_refinement_candidates_stay_inside_initial_set(&refinement);
        assert!(refinement.constraints().iter().any(|constraint| matches!(
            constraint,
            LocalEndomorphismConstraint::NodeLevel {
                provenance: ConstraintSource::NodeReport {
                    node_id: IsogenyGraphNodeId(0),
                },
                ..
            }
        )));
        assert_only_node_or_unambiguous_edge_constraints(&refinement);
    }

    #[test]
    fn refinement_rejects_missing_nodes() {
        let (_, report) = graph_report(0);

        let error = report
            .refine_candidates_for_node(
                IsogenyGraphNodeId(99),
                CandidateRefinementStrategy::NodeLocalLevelsOnly,
            )
            .expect_err("missing node should be rejected");

        assert_eq!(
            error,
            CandidateRefinementError::NodeNotFound {
                node_id: IsogenyGraphNodeId(99)
            }
        );
    }
}
