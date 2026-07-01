//! CM-inspired Frobenius trace candidates.
//!
//! This module contains the first small bridge from complex multiplication
//! heuristics to the Frobenius layer. For a curve with CM by a negative
//! discriminant `D`, the ordinary prime-field trace often satisfies
//! `4p = t² + |D|v²`. The helpers here only compute the arithmetic candidates
//! for `|t|`; they do not certify that a concrete curve has CM by `D`, and they
//! do not determine the sign of `t`.

mod trace_candidates;

#[cfg(test)]
mod tests;

pub use trace_candidates::{CmTraceCandidate, CmTraceCandidateError, cm_absolute_trace_candidates};
