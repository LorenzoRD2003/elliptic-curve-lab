//! Rational torsion scaffolding for short-Weierstrass curves over `Q`.
//!
//! This module will host the staged `E(Q)_tors` workflow:
//!
//! 1. transport a rational short-Weierstrass model to an integral model;
//! 2. enumerate integral candidates through Lutz-Nagell;
//! 3. verify candidates against Mazur's theorem;
//! 4. later compare that route with the reduction-mod-`p` and Hensel strategy
//!    from the problem-set exercise.
//!
//! The first milestone intentionally exposes only the shared value objects and
//! fixtures that later implementation passes will fill in.

mod classification;
mod enumeration;
mod error;
mod integral_model;
mod mazur;
mod report;
mod verification;

#[cfg(test)]
mod tests;

pub(crate) use classification::{RationalTorsionGroup, RationalTorsionGroupShape};
pub(crate) use error::RationalTorsionError;
pub(crate) use report::RationalTorsionReport;
