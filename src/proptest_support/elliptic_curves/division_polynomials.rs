use crate::fields::traits::*;
use proptest::prelude::*;
use std::fmt;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::division_polynomials::DivisionPolynomialForm;
use crate::fields::traits::EnumerableFiniteField;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Division-polynomial fixture derived from one short-Weierstrass curve.
#[derive(Clone)]
pub struct DivisionPolynomialCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub index: usize,
    pub polynomial: DivisionPolynomialForm<F>,
}

impl<F: Field> fmt::Debug for DivisionPolynomialCase<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DivisionPolynomialCase")
            .field("index", &self.index)
            .finish()
    }
}

/// Returns a division-polynomial case for small supported indices.
pub fn arb_division_polynomial_case<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<DivisionPolynomialCase<F>>
where
    F: EnumerableFiniteField + 'static,
    F::Elem: 'static,
{
    let max_index = config.max_division_index.max(1);

    arb_nonsingular_curve::<F>(config)
        .prop_flat_map(move |curve| {
            (Just(curve), 1usize..=max_index).prop_filter_map(
                "division polynomial should be supported at the sampled index",
                |(curve, index)| {
                    curve
                        .division_polynomial(index)
                        .ok()
                        .map(|polynomial| DivisionPolynomialCase {
                            curve,
                            index,
                            polynomial,
                        })
                },
            )
        })
        .boxed()
}

pub(crate) fn touch_division_polynomial_case_fields() {
    let _ = |case: DivisionPolynomialCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.index;
        let _ = case.polynomial;
    };
}
