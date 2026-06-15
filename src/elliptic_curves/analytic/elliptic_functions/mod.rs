//! Truncated Weierstrass elliptic functions attached to one lattice.
//!
//! For a lattice `Λ = ℤω₁ + ℤω₂`, the classical Weierstrass functions are
//! defined by lattice sums such as
//! `℘(z) = 1/z² + Σ_{ω ∈ Λ \\ {0}} (1/(z-ω)² - 1/ω²)`.
//! This module implements the first educational approximations to that story
//! by replacing the infinite lattice with a validated finite square box in
//! `ℤ²`.
//!
//! The current split is:
//!
//! - `truncation.rs` validates the finite lattice-sum policy.
//! - `traits.rs` and `evaluator.rs` keep the shared summation pattern local.
//! - `p.rs` evaluates truncated `℘`.
//! - `p_derivative.rs` evaluates truncated `℘′`.
//!
//! These approximations are meant to support later uniformization and
//! differential-equation checks, so the module emphasizes explicit reports and
//! validated truncation choices over raw numerical convenience.
mod evaluator;
mod p;
mod p_derivative;
pub(crate) mod traits;
mod truncation;

#[cfg(test)]
mod tests;

pub use p::WeierstrassPApprox;
pub use p_derivative::WeierstrassPDerivativeApprox;
pub use truncation::EllipticFunctionTruncation;
