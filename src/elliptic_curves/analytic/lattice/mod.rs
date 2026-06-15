//! Validated lattice and torus-side coordinates for analytic elliptic curves.
//!
//! The analytic story begins with a rank-two lattice
//! `Λ = ℤω₁ + ℤω₂ ⊂ ℂ` and the quotient torus `ℂ / Λ`.
//! This module owns the concrete value objects that make that starting point
//! explicit in code.
//!
//! The internal split is:
//!
//! - `types.rs` stores the main public lattice and torus values.
//! - `basis.rs` validates ordered period bases and derived lattice metadata.
//! - `coordinates.rs` and `points.rs` organize representatives in the
//!   fundamental parallelogram and in `ℂ / Λ`.
//! - `torus.rs` handles torus-side reduction/comparison helpers.
//! - `truncation.rs` validates finite square-box truncations for lattice sums.
//! - `context.rs` provides a small crate-private trait for reports that carry
//!   the ambient pair `(τ, Λ)`.
//!
//! This module is the ambient owner for later analytic routines such as
//! lattice invariants, truncated `℘` evaluation, modular comparisons, and
//! forward or inverse uniformization.
mod basis;
mod context;
mod coordinates;
mod points;
mod torus;
mod truncation;
mod types;

#[cfg(test)]
mod tests;

pub use truncation::LatticeSumTruncation;
pub use types::{
    ComplexLattice, ComplexModuloLatticeComparison, ComplexTorusPoint,
    FundamentalParallelogramCoordinate, LatticeIndexPoint,
};

pub(crate) use context::HasAnalyticLatticeContext;
