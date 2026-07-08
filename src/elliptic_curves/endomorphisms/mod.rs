//! Endomorphism-side arithmetic and candidate-set scaffolding.
//!
//! The current arithmetic surface is centered on imaginary quadratic orders.
//! Its radicand-based entrypoints intentionally support only negative integer
//! radicands `m < 0`, exposing the derived fundamental discriminant `D_K` and
//! maximal order `O_K` rather than a larger public quadratic-field API.

pub mod candidate_sets;
pub mod quadratic_ideals;
pub mod quadratic_orders;

#[cfg(test)]
mod tests;

pub use candidate_sets::{EndomorphismRingCandidateSet, EndomorphismRingReport};
