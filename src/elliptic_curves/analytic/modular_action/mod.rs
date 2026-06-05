//! Modular actions on the upper half-plane.
//!
//! The point of this submodule is not just that matrices in `SL_2(ℤ)` can
//! *change* a parameter `τ`, but that they change it in a very structured way:
//! different parameters related by the modular action describe the same
//! underlying complex torus, only written with a different lattice basis.
//!
//! If `Λ_τ = ℤ + ℤτ`, then for
//! `γ = [[a, b], [c, d]] ∈ SL_2(ℤ)` one has
//! `γ(τ) = (aτ + b) / (cτ + d)` and
//! `Λ_{γ(τ)} = (1 / (cτ + d)) Λ_τ`.
//!
//! So the transformed lattice is not literally equal to the original one, but
//! it differs only by multiplication by a nonzero complex scalar. That means
//! `τ` and `γ(τ)` encode the same complex torus `ℂ / Λ` up to isomorphism.
//! This is why modular invariants such as `j` are expected to agree on both
//! sides: they depend on the geometric object represented by `τ`.
//!
//! In that language, the upper half-plane `ℍ` is a parameter space
//! with redundancy, and the actual moduli space of complex tori is modeled by
//! the quotient `SL_2(ℤ) \ ℍ`.
//! A single point `τ ∈ H` is one coordinate choice, while its whole
//! `SL_2(ℤ)`-orbit represents the same geometric torus.

mod invariance;
mod matrix;

pub use invariance::{ModularInvarianceReport, verify_j_modular_invariance};
pub use matrix::ModularMatrix;

#[cfg(test)]
mod tests;
