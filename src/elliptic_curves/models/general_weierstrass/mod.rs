//! General Weierstrass curves introduced in staged milestones.
//!
//! The long-term target is the affine model
//!
//! `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`.
//!
//! The current milestone intentionally exposes only the model descriptor and
//! its coefficients so the wider crate can start acknowledging a second curve
//! family without disturbing the mature short-Weierstrass stack.

mod type_definition;

pub use type_definition::GeneralWeierstrassCurve;
