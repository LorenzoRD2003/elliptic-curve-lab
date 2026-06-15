//! Truncated lattice Eisenstein sums attached to a complex lattice.
//!
//! This namespace keeps the approximation reports public while making
//! `ComplexLattice` the owner of the actual evaluation routes.

mod api;
mod types;

#[cfg(test)]
mod tests;

pub use types::{EisensteinSumApprox, TruncationConvergenceReport};
