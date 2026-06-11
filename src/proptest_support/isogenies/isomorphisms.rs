use proptest::prelude::*;

use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassIsomorphism};
use crate::fields::Fp;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::arb_curve_and_point;
use crate::proptest_support::fields::arb_nonzero_fp_elem;

/// Small isomorphism fixture over a short-Weierstrass curve.
#[derive(Clone, Debug)]
pub struct ShortWeierstrassIsomorphismCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub isomorphism: ShortWeierstrassIsomorphism<Fp<P>>,
    pub sample_point: crate::elliptic_curves::AffinePoint<Fp<P>>,
}

/// Returns a short-Weierstrass scaling isomorphism together with one point in
/// its domain.
pub fn arb_short_weierstrass_isomorphism_case<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<ShortWeierstrassIsomorphismCase<P>> {
    (arb_curve_and_point::<P>(config), arb_nonzero_fp_elem::<P>())
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
    let _ = |case: ShortWeierstrassIsomorphismCase<17>| {
        let _ = case.curve;
        let _ = case.isomorphism;
        let _ = case.sample_point;
    };
}
