//! Legendre normalization of the recovered Weierstrass cubic.
//!
//! Starting from three roots `e₁, e₂, e₃` of
//! `4x³ - g₂x - g₃ = 4(x-e₁)(x-e₂)(x-e₃)`,
//! one may transport the curve to the Legendre form
//! `Y² = X(X-1)(X-λ)`.
//! The choice is not unique: permuting the roots changes `λ` inside its
//! six-element `S₃` orbit.
//!
//! This module separates that story into:
//!
//! - `parameter.rs` for the validated Legendre parameter itself.
//! - `orbit.rs` for the full `S₃` orbit data.
//! - `conditioning.rs` for near-singularity diagnostics.
//! - `candidate.rs` for the internal deterministic candidate comparison.
//! - `reduction.rs` for the chosen affine transport and scale factors.
//! - `report.rs` for the structured explanation of one selected reduction.
//!
//! In other words, this module owns the passage from raw cubic roots to one
//! deterministic normalized Legendre model together with the diagnostics that
//! explain why that representative was chosen.
mod candidate;
mod conditioning;
mod orbit;
mod parameter;
mod reduction;
mod report;

pub use conditioning::LegendreParameterConditioning;
pub use orbit::{LegendreOrbitElement, LegendreOrbitElementKind, LegendreParameterOrbit};
pub use parameter::LegendreParameter;
pub use reduction::LegendreReduction;
pub use report::LegendreReductionReport;

pub(crate) use candidate::LegendreCandidate;
pub(crate) use reduction::LEGENDRE_PERMUTATIONS;
