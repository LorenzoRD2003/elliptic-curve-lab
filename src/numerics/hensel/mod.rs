//! Exact Hensel lifting helpers for integer-polynomial congruences.
//!
//! The general route lifts a simple root of a polynomial from modulo `p^k` to
//! modulo `p^(k + 1)` using the correction `x_{k+1} = x_k + t p^k`.
//! Specialized odd-prime square-root helpers solve `x² = a` by either
//! repeated adjacent lifts or by the fast Newton-Hensel route
//! `p → p² → p⁴ → ... → p^e`. The `2^e` square-root helper uses a
//! separate bit-lifting route because `2x` is never a unit modulo `2`.
//! The general integer-modulus square-root helper factors `m = Π pᵢ^eᵢ`,
//! solves each prime-power component, and recombines the local roots with CRT.
//!
//! All current entry points are crate-internal. The general and odd-prime
//! square-root routes assume the simple-root condition, such as
//! `f'(x) != 0 mod p` or `2x != 0 mod p`; the `2^e` route handles its singular
//! lifting separately.

#![cfg_attr(not(test), allow(dead_code, unused_imports))]

mod api;
mod error;
mod integer_roots;
mod polynomial;
mod square_root;
mod step;
mod trace;

#[cfg(test)]
mod tests;

pub(crate) use api::{hensel_lift_simple_root, hensel_lift_simple_root_step};
pub(crate) use error::HenselLiftError;
pub(crate) use integer_roots::{
    HenselIntegerRootSearchConfig, HenselIntegerRootSearchReport, HenselIntegerRootTrace,
    find_integer_roots_by_hensel, hensel_lift_integer_root,
};
pub(crate) use square_root::{
    hensel_lift_square_root, hensel_lift_square_root_fast, sqrt_mod_m, sqrt_mod_odd_prime_power,
    sqrt_mod_odd_prime_power_divisible_radicand, sqrt_mod_two_power,
};

#[allow(unused_imports)]
pub(crate) use step::{HenselLiftStep, HenselSquareRootFastStep};
#[allow(unused_imports)]
pub(crate) use trace::{HenselLiftTrace, HenselSquareRootFastTrace};
