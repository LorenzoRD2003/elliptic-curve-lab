use core::fmt::Debug;

use proptest::prelude::*;

use crate::fields::{rational_function_field::RationalFunction, traits::EnumerableFiniteField};
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::{arb_dense_polynomial, arb_nonzero_dense_polynomial};

/// Returns a strategy for canonical rational functions over a small
/// enumerable finite field.
///
/// The generator samples dense numerators and non-zero dense denominators, then
/// relies on [`RationalFunction::new`] to reduce the presentation to the
/// repository's canonical rational-function form.
pub fn arb_rational_function<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<RationalFunction<F>>
where
    F: EnumerableFiniteField + Debug + 'static,
{
    (
        arb_dense_polynomial::<F>(config),
        arb_nonzero_dense_polynomial::<F>(config),
    )
        .prop_map(|(numerator, denominator)| {
            RationalFunction::<F>::new(numerator, denominator)
                .expect("non-zero denominator should define a rational function")
        })
        .boxed()
}
