//! Exact Hensel lifting helpers for integer-polynomial congruences.
//!
//! The general route lifts a simple root of a polynomial from modulo `p^k` to
//! modulo `p^(k + 1)` using the correction `x_{k+1} = x_k + t p^k`.
//! Specialized square-root helpers solve `x^2 = a` by either repeated adjacent
//! lifts or by the fast Newton-Hensel route
//! `p -> p^2 -> p^4 -> ... -> p^e`.
//!
//! All current entry points are crate-internal and assume the simple-root
//! condition, such as `f'(x) != 0 mod p` or `2x != 0 mod p` for square roots.
//! Singular Hensel lifting has different branching behavior and should grow as
//! a separate surface when the repo needs it.

#![cfg_attr(not(test), allow(dead_code, unused_imports))]

mod api;
mod error;
mod polynomial;
mod square_root;
mod step;
mod trace;

#[cfg(test)]
mod tests;

pub(crate) use api::{hensel_lift_simple_root, hensel_lift_simple_root_step};
pub(crate) use error::HenselLiftError;
pub(crate) use square_root::{hensel_lift_square_root, hensel_lift_square_root_fast};
#[allow(unused_imports)]
pub(crate) use step::{HenselLiftStep, HenselSquareRootFastStep};
#[allow(unused_imports)]
pub(crate) use trace::{HenselLiftTrace, HenselSquareRootFastTrace};
