use proptest::prelude::*;

use crate::elliptic_curves::frobenius::{FrobeniusDiscriminant, FrobeniusTrace};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::FrobeniusTraceCurveModel;
use crate::fields::Fp;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Frobenius-side case packaged from one finite short-Weierstrass curve.
#[derive(Clone, Debug)]
pub struct FrobeniusCurveCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub trace: FrobeniusTrace,
    pub discriminant: FrobeniusDiscriminant,
}

/// Returns a finite-field curve together with its Frobenius trace package.
pub fn arb_frobenius_curve_case<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<FrobeniusCurveCase<P>> {
    arb_nonsingular_curve::<P>(config)
        .prop_map(|curve| {
            let trace = curve
                .frobenius_trace()
                .expect("enumerable finite short-Weierstrass curves should admit a trace");
            let discriminant = trace.discriminant();

            FrobeniusCurveCase {
                curve,
                trace,
                discriminant,
            }
        })
        .boxed()
}

pub(crate) fn touch_frobenius_case_fields() {
    let _ = |case: FrobeniusCurveCase<17>| {
        let _ = case.curve;
        let _ = case.trace;
        let _ = case.discriminant;
    };
}
