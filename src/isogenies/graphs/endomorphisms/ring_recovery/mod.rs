//! Endomorphism-ring level recovery from volcano floor distances.
//!
//! The local report recovers one exponent `v_ℓ(u)`. The global report only
//! assembles those already-certified local exponents; it does not coordinate
//! construction or search across several `ℓ`-graphs.

mod error;
mod global;
mod local;

#[cfg(test)]
mod tests;

pub use error::EndomorphismRingLevelRecoveryError;
pub use global::EndomorphismRingLevelRecoveryReport;
pub use local::LocalEndomorphismRingLevelReport;
