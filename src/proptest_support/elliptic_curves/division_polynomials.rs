use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::division_polynomials::{DivisionPolynomial, division_polynomial};
use crate::fields::Fp;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Division-polynomial fixture derived from one short-Weierstrass curve.
#[derive(Clone, Debug)]
pub struct DivisionPolynomialCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub index: usize,
    pub polynomial: DivisionPolynomial<Fp<P>>,
}

/// Returns a division-polynomial case for small supported indices.
pub fn arb_division_polynomial_case<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<DivisionPolynomialCase<P>> {
    let max_index = config.max_division_index.max(1);

    arb_nonsingular_curve::<P>(config)
        .prop_flat_map(move |curve| {
            (Just(curve), 1usize..=max_index).prop_filter_map(
                "division polynomial should be supported at the sampled index",
                |(curve, index)| {
                    division_polynomial(&curve, index).ok().map(|polynomial| {
                        DivisionPolynomialCase {
                            curve,
                            index,
                            polynomial,
                        }
                    })
                },
            )
        })
        .boxed()
}

pub(crate) fn touch_division_polynomial_case_fields() {
    let _ = |case: DivisionPolynomialCase<17>| {
        let _ = case.curve;
        let _ = case.index;
        let _ = case.polynomial;
    };
}
