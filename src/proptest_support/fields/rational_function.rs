use core::fmt::Debug;

use proptest::prelude::*;

use crate::fields::{EnumerableFiniteField, RationalFunction};
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::arb_dense_polynomial;

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
    let denominator_config = PolynomialStrategyConfig {
        require_nonzero_leading_coefficient: true,
        ..config
    };

    (
        arb_dense_polynomial::<F>(config),
        arb_dense_polynomial::<F>(denominator_config).prop_filter(
            "rational-function denominator should be non-zero",
            |polynomial| !polynomial.is_zero(),
        ),
    )
        .prop_map(|(numerator, denominator)| {
            RationalFunction::<F>::new(numerator, denominator)
                .expect("non-zero denominator should define a rational function")
        })
        .boxed()
}
