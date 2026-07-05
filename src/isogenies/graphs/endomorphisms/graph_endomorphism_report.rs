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
    /// In the current staged implementation, [`CandidateRefinementStrategy::Conservative`]
    /// and [`CandidateRefinementStrategy::NodeLocalLevelsOnly`] both use only
    /// the node-local conductor-level constraint `v_ℓ(f) ∈ L`. Incident edge
    /// evidence is reserved for the next refinement phase.
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
            CandidateRefinementStrategy::Conservative
            | CandidateRefinementStrategy::NodeLocalLevelsOnly => {
                EndomorphismCandidateRefinement::from_node_local_levels(node_report, self.prime())
            }
            CandidateRefinementStrategy::IncidentUnambiguousEdges => {
                Err(CandidateRefinementError::StrategyNotImplemented { strategy })
            }
        }
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
            IsogenyEdgeEndomorphismTentativeRelation,
            refinement::{
                CandidateRefinementError, CandidateRefinementStrategy, ConstraintSource,
                LocalEndomorphismConstraint, RefinementConfidence,
            },
        },
    };

    type F41 = crate::fields::Fp41;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn graph_report_collects_node_and_edge_endomorphism_data() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

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
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

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
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

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
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");
        let node_report = report
            .node_report(IsogenyGraphNodeId(0))
            .expect("root node report should exist");

        let refinement = report
            .refine_candidates_for_node(
                IsogenyGraphNodeId(0),
                CandidateRefinementStrategy::NodeLocalLevelsOnly,
            )
            .expect("node-local refinement should build");

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
        assert_eq!(ell, &BigUint::from(2u8));
        assert_eq!(
            allowed_levels,
            &node_report.possible_levels().into_iter().collect()
        );
        assert_eq!(
            provenance,
            &ConstraintSource::NodeReport {
                node_id: IsogenyGraphNodeId(0)
            }
        );
    }

    #[test]
    fn conservative_refinement_currently_uses_node_local_evidence() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

        let conservative = report
            .refine_candidates_for_node(
                IsogenyGraphNodeId(0),
                CandidateRefinementStrategy::Conservative,
            )
            .expect("conservative refinement should build");
        let node_local = report
            .refine_candidates_for_node(
                IsogenyGraphNodeId(0),
                CandidateRefinementStrategy::NodeLocalLevelsOnly,
            )
            .expect("node-local refinement should build");

        assert_eq!(conservative, node_local);
    }

    #[test]
    fn incident_edge_refinement_is_explicitly_deferred() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

        let error = report
            .refine_candidates_for_node(
                IsogenyGraphNodeId(0),
                CandidateRefinementStrategy::IncidentUnambiguousEdges,
            )
            .expect_err("incident-edge refinement should be phase-3 work");

        assert_eq!(
            error,
            CandidateRefinementError::StrategyNotImplemented {
                strategy: CandidateRefinementStrategy::IncidentUnambiguousEdges
            }
        );
    }

    #[test]
    fn refinement_rejects_missing_nodes() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("graph endomorphism report should build");

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
