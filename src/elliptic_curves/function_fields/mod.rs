//! Function fields of short-Weierstrass elliptic curves.
//!
//! This module currently implements the function field of a validated short
//! Weierstrass curve `E : y^2 = x^3 + ax + b` over a base field `F` of
//! characteristic different from `2` and `3`.
//!
//! The implementation uses the degree-two presentation `F(E) = F(x) ⊕ y F(x)`
//! as a vector space over the rational function field `F(x)`. Concretely, an
//! element is stored as a pair of rational functions `(A, B)` representing
//! `A(x) + y B(x)`.
//!
//! The module also exposes a small point representation for curve points with
//! coordinates in `F(E)`. Its first purpose is pedagogical: it models the
//! generic point `(x, y)` and the rational group-law formulas applied to that
//! point, which is the natural route to pullback formulas such as `[n]^*(x)`
//! and `[n]^*(y)`.
//!
//! Multiplication is reduced using the short-Weierstrass relation
//! `y^2 = f(x) = x^3 + ax + b`,
//! so `(A, B) * (C, D) = (AC + fBD, AD + BC)`.
//!
//! The conjugation involution sends `y` to `-y`, hence `conj(A, B) = (A, -B)`,
//! and the norm is `N(A, B) = A^2 - fB^2`. When `N(A, B)` is non-zero, the
//! inverse is computed by the classical conjugate-over-norm formula
//! `(A, B)^(-1) = (A / N(A, B), -B / N(A, B))`.

mod field;
mod point;
mod value;

#[cfg(test)]
mod tests;

pub use field::ShortWeierstrassFunctionField;
pub use point::ShortWeierstrassFunctionFieldPoint;
pub use value::ShortWeierstrassFunction;
