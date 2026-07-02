use crate::fields::traits::*;
use proptest::prelude::*;

use crate::fields::traits::EnumerableFiniteField;
use crate::polynomials::DensePolynomial;
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::shared::{arb_field_elem, arb_nonzero_field_elem};

/// Returns a dense polynomial over a small enumerable finite field.
pub fn arb_dense_polynomial<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<DensePolynomial<F>>
where
    F: EnumerableFiniteField + 'static,
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

/// Returns a non-zero dense polynomial over a small enumerable finite field.
///
/// The sample is constructed with a non-zero trailing coefficient instead of
/// relying on rejection, which makes it a better fit for denominators and
/// other contexts that require a genuine leading term.
pub fn arb_nonzero_dense_polynomial<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<DensePolynomial<F>>
where
    F: EnumerableFiniteField + 'static,
{
    if config.max_len == 0 {
        return Just(DensePolynomial::<F>::new(Vec::new())).boxed();
    }

    (0usize..config.max_len)
        .prop_flat_map(move |prefix_len| {
            (
                prop::collection::vec(arb_field_elem::<F>(), prefix_len),
                arb_nonzero_field_elem::<F>(),
            )
                .prop_map(|(mut coefficients, leading)| {
                    coefficients.push(leading);
                    DensePolynomial::<F>::new(coefficients)
                })
        })
        .boxed()
}
