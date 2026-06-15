//! Analytic invariants attached to lattices and upper-half-plane parameters.
//!
//! This namespace keeps the value objects public while attaching the main
//! computations to their natural owners: lattices, upper-half-plane points,
//! and the invariant bundle itself.

mod api;
mod value;

#[cfg(test)]
mod special_values;
#[cfg(test)]
mod tests;

pub use value::AnalyticInvariants;

#[cfg(test)]
pub(crate) use special_values::{ComplexAnalyticCurveLabReport, SpecialJKind, SpecialTauKind};
