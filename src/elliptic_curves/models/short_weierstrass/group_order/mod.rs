//! Curve-side group-order entry points for short-Weierstrass models.
//!
//! The public report and strategy types currently live under
//! `elliptic_curves::frobenius::group_order`, because they package Frobenius
//! data such as `#E(F_q)`, `t`, and `H(q)` in a family-independent way.
//!
//! This `short_weierstrass::group_order` module is therefore the curve-side
//! owner of the executable algorithms and dispatch:
//! - deterministic routes in [`api`]
//! - the quadratic-character counting route in [`quadratic_character`]
//! - the automatic Schoof route in [`api`]
//! - Mestre's prime-field route in [`mestre`]

mod api;
mod mestre;
mod quadratic_character;

#[cfg(test)]
mod tests;

pub use crate::elliptic_curves::frobenius::group_order::{
    FiniteFieldGroupOrderStrategy, GroupOrderReport, GroupOrderRoute, MestreConfig,
    MestreGroupOrderReport, SchoofGroupOrderSummary, SmallFieldGroupOrderStrategy,
    SmallFieldSampledGroupOrderStrategy,
};
