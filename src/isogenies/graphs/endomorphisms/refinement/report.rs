use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::ImaginaryQuadraticOrder,
};
use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::refinement::{
        CandidateElimination, CandidateRefinementStrategy, LocalEndomorphismConstraint,
        RefinementConfidence,
    },
};

/// One node's candidate-refinement run.
///
/// A refinement run starts from the Frobenius-compatible set
/// `C₀ = {O_f : ℤ[π] ⊆ O_f ⊆ O_K}` and records the survivor set
/// `C_surv ⊆ C₀` after applying local graph evidence. Even when
/// `C_surv` has one element, the result means “uniquely compatible with the
/// observed evidence”, not “certified as `End(E)`”.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndomorphismCandidateRefinement {
    node_id: IsogenyGraphNodeId,
    initial_candidates: EndomorphismRingCandidateSet,
    surviving_candidates: Vec<ImaginaryQuadraticOrder>,
    eliminated_candidates: Vec<CandidateElimination>,
    constraints: Vec<LocalEndomorphismConstraint>,
    confidence: RefinementConfidence,
}

impl EndomorphismCandidateRefinement {
    pub(crate) fn new(
        node_id: IsogenyGraphNodeId,
        initial_candidates: EndomorphismRingCandidateSet,
        surviving_candidates: Vec<ImaginaryQuadraticOrder>,
        eliminated_candidates: Vec<CandidateElimination>,
        constraints: Vec<LocalEndomorphismConstraint>,
        confidence: RefinementConfidence,
    ) -> Self {
        Self {
            node_id,
            initial_candidates,
            surviving_candidates,
            eliminated_candidates,
            constraints,
            confidence,
        }
    }

    /// Returns the node refined by this run.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the initial Frobenius-compatible candidate set `C₀`.
    pub fn initial_candidates(&self) -> &EndomorphismRingCandidateSet {
        &self.initial_candidates
    }

    /// Returns the candidate orders that survived the applied evidence.
    pub fn surviving_candidates(&self) -> &[ImaginaryQuadraticOrder] {
        &self.surviving_candidates
    }

    /// Returns the candidates removed during the refinement run.
    pub fn eliminated_candidates(&self) -> &[CandidateElimination] {
        &self.eliminated_candidates
    }

    /// Returns the local constraints applied by this refinement run.
    pub fn constraints(&self) -> &[LocalEndomorphismConstraint] {
        &self.constraints
    }

    /// Returns the kind of evidence supporting the survivor set.
    pub fn confidence(&self) -> RefinementConfidence {
        self.confidence
    }

    /// Returns the unique surviving candidate when the evidence leaves
    /// exactly one `O_f`.
    pub fn unique_survivor(&self) -> Option<&ImaginaryQuadraticOrder> {
        if self.surviving_candidates.len() == 1 {
            self.surviving_candidates.first()
        } else {
            None
        }
    }

    /// Returns whether the observed evidence leaves exactly one surviving
    /// candidate.
    pub fn is_unique(&self) -> bool {
        self.unique_survivor().is_some()
    }
}

/// Aggregate refinement report for every node in one graph endomorphism report.
///
/// This is the graph-level counterpart of
/// [`EndomorphismCandidateRefinement`]. It packages one refinement run per
/// node for a fixed local prime `ℓ` and strategy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphCandidateRefinementReport {
    prime: BigUint,
    strategy: CandidateRefinementStrategy,
    node_refinements: Vec<EndomorphismCandidateRefinement>,
}

impl IsogenyGraphCandidateRefinementReport {
    /// Returns the chosen local prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the strategy used to construct the node refinements.
    pub fn strategy(&self) -> CandidateRefinementStrategy {
        self.strategy
    }

    /// Returns the node refinements in dense node-id order.
    pub fn node_refinements(&self) -> &[EndomorphismCandidateRefinement] {
        &self.node_refinements
    }

    /// Returns the refinement for the requested node when present.
    pub fn refinement_for_node(
        &self,
        node_id: IsogenyGraphNodeId,
    ) -> Option<&EndomorphismCandidateRefinement> {
        self.node_refinements
            .get(node_id.0)
            .filter(|refinement| refinement.node_id == node_id)
    }
}
