use core::fmt::Debug;

use proptest::prelude::*;

use crate::fields::EnumerableFiniteField;
use crate::polynomials::DensePolynomial;
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::shared::arb_field_elem;

/// Returns a dense polynomial over a small enumerable finite field.
pub fn arb_dense_polynomial<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<DensePolynomial<F>>
where
    F: EnumerableFiniteField + Debug + 'static,
{
    prop::collection::vec(arb_field_elem::<F>(), 0..=config.max_len)
        .prop_map(DensePolynomial::<F>::new)
        .prop_filter(
            "dense polynomial should keep a non-zero leading coefficient when requested",
            move |polynomial| {
                !config.require_nonzero_leading_coefficient
                    || polynomial
                        .coefficients()
                        .last()
                        .is_none_or(|coefficient| !F::is_zero(coefficient))
            },
        )
        .boxed()
}
