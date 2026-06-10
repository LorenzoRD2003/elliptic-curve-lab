//! Polynomial-oriented strategies organized by representation.

pub mod dense;
pub mod multivariate;
pub mod shared;
pub mod sparse;

pub use dense::arb_dense_polynomial;
pub use multivariate::arb_multivariate_polynomial;
pub use sparse::arb_sparse_polynomial;

pub(crate) fn touch_polynomial_inventory() {
    let config = crate::proptest_support::config::PolynomialStrategyConfig::default();
    let _ = arb_dense_polynomial::<crate::fields::Fp<17>>(config);
    let _ = arb_sparse_polynomial::<crate::fields::Fp<17>>(config);
    let _ = arb_multivariate_polynomial::<crate::fields::Fp<17>>(config);
}
