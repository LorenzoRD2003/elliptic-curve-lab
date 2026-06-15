//! Short-Weierstrass-specific isogeny constructions and reports.
//!
//! This namespace owns the concrete isogeny stories that are currently
//! implemented only for short-Weierstrass curves:
//!
//! - classical Vélu constructions from explicit finite kernels;
//! - pullbacks on short-Weierstrass function fields;
//! - absolute/relative Frobenius isogenies and Verschiebung certificates;
//! - short-Weierstrass duality builders and reports.

mod curve_api;
pub mod frobenius;
pub mod function_field_maps;
pub mod velu;

pub use crate::isogenies::velu::VeluIsogeny;
pub use velu::DualVeluIsogeny;
