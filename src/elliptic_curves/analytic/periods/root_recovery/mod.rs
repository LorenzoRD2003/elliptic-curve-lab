//! Recovery of the Weierstrass cubic roots from analytic invariants.
//!
//! For an analytic curve `y² = 4x³ - g₂x - g₃`, one wants approximate roots
//! `e₁, e₂, e₃` of the cubic factor
//! `4(x-e₁)(x-e₂)(x-e₃)`.
//! The current implementation uses a hybrid route:
//!
//! - a Cardano-style rough recovery away from unstable regimes,
//! - a special near-pure-cubic fallback when the depressed cubic is
//!   numerically close to `x³ + q`,
//! - Newton polishing on each recovered root,
//! - explicit reconstruction checks against `g₂`, `g₃`, and the vanishing
//!   `x²` coefficient.
//!
//! The files are split accordingly:
//!
//! - `api.rs` exposes the curve-level recovery methods.
//! - `internal.rs` contains the hybrid recovery logic and polishing helpers.
//! - `report.rs` stores the structured recovery report and Cardano-specific
//!   diagnostics.
//!
//! This is an approximate educational routine, so branch-choice and recovery
//! diagnostics remain visible rather than being collapsed to a single bare
//! root triple.
mod api;
mod internal;
mod report;

pub use report::{CardanoRootRecoveryDiagnostics, CubicRootRecoveryReport};
