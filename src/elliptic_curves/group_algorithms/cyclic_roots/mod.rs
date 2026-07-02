//! Prime-degree root extraction in externally certified finite cyclic groups.
//!
//! The public surface is intentionally the curve method
//! [`CyclicGroupPrimeRootCurveModel::cyclic_group_prime_root`] plus read-only
//! reports. The executable engine stays internal so callers cannot bypass the
//! curve-side validation story.

#![cfg_attr(not(test), allow(dead_code))]

mod algorithm;
mod bezout;
mod curve_model;
mod error;
mod input;
mod outcome;
mod report;
mod step;
mod trace;

#[cfg(test)]
mod tests;

pub use bezout::CyclicPrimeRootBezout;
pub use curve_model::CyclicGroupPrimeRootCurveModel;
pub use error::CyclicPrimeRootError;
pub use input::{CyclicPrimeRootInput, CyclicPrimeRootInputError};
pub use outcome::CyclicPrimeRootOutcome;
pub use report::CyclicPrimeRootReport;
pub use step::CyclicPrimeRootStep;
pub use trace::CyclicPrimeRootTrace;
