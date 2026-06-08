use std::hash::Hash;

use crate::elliptic_curves::endomorphisms::{
    EndomorphismRingCandidateSet, VolcanoEndomorphismLevelCandidate,
};
use crate::elliptic_curves::frobenius::FrobeniusTraceCurveModel;
use crate::fields::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyEdgeEndomorphismReport, IsogenyGraph, IsogenyGraphEdgeId,
    IsogenyGraphError, IsogenyGraphNodeId,
};
use num_bigint::BigUint;

/// Endomorphism-side report for one stored graph node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismNodeReport {
    node_id: IsogenyGraphNodeId,
    candidate_set: EndomorphismRingCandidateSet,
    local_levels: Vec<VolcanoEndomorphismLevelCandidate>,
    possible_levels: Vec<u32>,
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
    /// Returns the node identifier.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the Frobenius-compatible candidate orders for this node.
    pub fn candidate_set(&self) -> &EndomorphismRingCandidateSet {
        &self.candidate_set
    }

    /// Returns the `ℓ`-local level candidates attached to this node.
    pub fn local_levels(&self) -> &[VolcanoEndomorphismLevelCandidate] {
        &self.local_levels
    }

    /// Returns the distinct possible local levels for this node.
    pub fn possible_levels(&self) -> &[u32] {
        &self.possible_levels
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
    pub fn from_graph<C: GraphCurveModel + FrobeniusTraceCurveModel>(
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
                let possible_levels = distinct_levels(&local_levels);

                Ok(IsogenyGraphEndomorphismNodeReport {
                    node_id: *node_id,
                    candidate_set: candidate_set.clone(),
                    local_levels,
                    possible_levels,
                })
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

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::{
        IsogenyEdgeEndomorphismRelation, IsogenyGraphBuilder, IsogenyGraphEndomorphismReport,
        IsogenyGraphNodeId,
    };
    use num_bigint::BigUint;

    type F41 = Fp<41>;
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

        let report = IsogenyGraphEndomorphismReport::from_graph(&graph, &BigUint::from(2u8))
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
        assert!(
            report
                .nodes()
                .iter()
                .all(|node| !node.local_levels().is_empty() && !node.possible_levels().is_empty())
        );
    }

    #[test]
    fn graph_report_exposes_dense_node_lookup() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(0)
            .build()
            .expect("depth-zero graph should build");

        let report = IsogenyGraphEndomorphismReport::from_graph(&graph, &BigUint::from(2u8))
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

        let report = IsogenyGraphEndomorphismReport::from_graph(&graph, &BigUint::from(2u8))
            .expect("graph endomorphism report should build");

        assert!(report.edges().iter().all(|edge| {
            matches!(
                edge.relation().relation(),
                IsogenyEdgeEndomorphismRelation::PossiblyHorizontal
                    | IsogenyEdgeEndomorphismRelation::PossiblyAscending
                    | IsogenyEdgeEndomorphismRelation::PossiblyDescending
                    | IsogenyEdgeEndomorphismRelation::Ambiguous
                    | IsogenyEdgeEndomorphismRelation::Unsupported
            )
        }));
    }
}
