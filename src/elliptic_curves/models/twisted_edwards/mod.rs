//! Twisted Edwards curves introduced in staged milestones.
//!
//! The current milestone supports:
//!
//! - validated curve descriptors for
//!   `a x^2 + y^2 = 1 + d x^2 y^2`
//! - classical invariants
//!
//! over fields of characteristic different from `2`.

mod invariants;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::TwistedEdwardsCurve;
