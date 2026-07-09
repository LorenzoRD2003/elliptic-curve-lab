//! Integral binary quadratic forms for the endomorphism-side class story.
//!
//! A form is stored by its canonical integral ternary coefficients
//! `(a, b, c)` and represents `ax² + bxy + cy²`.
//!
//! The multivariate-polynomial view is derived over `Q`, because the polynomial
//! engine is field-based while binary quadratic forms are integral objects.

mod form;

#[cfg(test)]
mod tests;

pub use form::BinaryQuadraticForm;
