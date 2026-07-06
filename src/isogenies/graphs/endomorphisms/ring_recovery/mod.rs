//! Endomorphism-ring level recovery from volcano floor distances.
//!
//! [`IsogenyGraph::recover_endomorphism_ring_at`](crate::isogenies::graphs::IsogenyGraph::recover_endomorphism_ring_at)
//! is the user-facing coordinator: it derives the Frobenius-compatible
//! candidates for one graph node, builds the required local `ℓ`-graphs, and
//! returns one assembled report. The local report still recovers one exponent
//! `v_ℓ(u)`, while the global report remains the validation/assembly surface
//! for already-certified local exponents.

mod error;
mod global;
mod local;
mod recover;

#[cfg(test)]
mod tests;

pub use error::EndomorphismRingLevelRecoveryError;
pub use global::EndomorphismRingLevelRecoveryReport;
pub use local::LocalEndomorphismRingLevelReport;
