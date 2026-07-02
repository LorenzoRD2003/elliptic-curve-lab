use crate::fields::traits::*;
use proptest::prelude::*;
use std::fmt;

use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::isomorphisms::ShortWeierstrassIsomorphism;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::arb_curve_and_point;
use crate::proptest_support::fields::arb_nonzero_fp_elem;

/// Small isomorphism fixture over a short-Weierstrass curve.
#[derive(Clone)]
pub struct ShortWeierstrassIsomorphismCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub isomorphism: ShortWeierstrassIsomorphism<F>,
    pub sample_point: crate::elliptic_curves::AffinePoint<F>,
}

impl<F: Field> fmt::Debug for ShortWeierstrassIsomorphismCase<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ShortWeierstrassIsomorphismCase")
            .finish()
    }
}

/// Returns a short-Weierstrass scaling isomorphism together with one point in
/// its domain.
pub fn arb_short_weierstrass_isomorphism_case<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<ShortWeierstrassIsomorphismCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
{
    (arb_curve_and_point::<F>(config), arb_nonzero_fp_elem::<F>())
        .prop_map(|((curve, sample_point), scaling_factor)| {
            let isomorphism = ShortWeierstrassIsomorphism::new(curve.clone(), scaling_factor)
                .expect("non-zero scaling factor should define a scaling isomorphism");
            ShortWeierstrassIsomorphismCase {
                curve,
                isomorphism,
                sample_point,
            }
        })
        .boxed()
}

pub(crate) fn touch_isomorphism_case_fields() {
    let _ = |case: ShortWeierstrassIsomorphismCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.isomorphism;
        let _ = case.sample_point;
    };
}
