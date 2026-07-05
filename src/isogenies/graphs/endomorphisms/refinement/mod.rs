//! Candidate-refinement surfaces for graph-side endomorphism evidence.
//!
//! This module records the vocabulary for moving from the Frobenius-compatible
//! candidate set `C₀ = {O_f: ℤ[π] ⊆ O_f ⊆ O_K}` to a survivor set
//! `C_surv ⊆ C₀` using local `ℓ`-isogeny-graph evidence.
//!
//! The types here are deliberately report-oriented. They do not certify
//! `End(E)`; they describe which candidates remain compatible with the
//! evidence observed by an [`IsogenyGraphEndomorphismReport`].
//!
//! [`IsogenyGraphEndomorphismReport`]: crate::isogenies::graphs::IsogenyGraphEndomorphismReport

mod confidence;
mod constraint;
mod elimination;
mod error;
mod report;
mod strategy;

pub use confidence::RefinementConfidence;
pub use constraint::{ConstraintSource, LocalEndomorphismConstraint};
pub use elimination::{
    CandidateElimination, CandidateEliminationReason, CandidateRefinementEdgeDirection,
};
pub use error::CandidateRefinementError;
pub use report::{EndomorphismCandidateRefinement, IsogenyGraphCandidateRefinementReport};
pub use strategy::CandidateRefinementStrategy;
