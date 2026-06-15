//! Forward uniformization from torus representatives to the analytic cubic.
//!
//! Once a lattice `Λ` is fixed, the classical map
//! `z mod Λ ↦ (℘(z), ℘′(z))`
//! sends torus classes to points on the analytic Weierstrass curve, with
//! lattice points mapping to the point at infinity.
//! This module contains the first explicit forward side of that story.
//!
//! The current split is:
//!
//! - `forward_map.rs` for the actual torus-to-curve map and its values.
//! - `differential_equation.rs` for checking the defining equation
//!   `℘′² = 4℘³ - g₂℘ - g₃`.
//! - `experiment.rs` for higher-level bundled reports used in examples and
//!   visualizations.
//!
//! This module is intentionally the forward companion to
//! `inverse_uniformization/`: one maps torus classes to the curve here, and
//! later tries to recover torus classes back from curve points there.
mod differential_equation;
mod experiment;
mod forward_map;

#[cfg(test)]
mod tests;

pub use differential_equation::{
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
};
pub use forward_map::{TorusToCurveMapResult, TorusToCurveValues};

#[cfg(test)]
pub(crate) use experiment::UniformizationExperimentReport;
