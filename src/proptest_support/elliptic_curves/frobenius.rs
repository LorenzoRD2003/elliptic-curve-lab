use crate::fields::traits::*;
use proptest::prelude::*;
use std::fmt;

use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::frobenius::{FrobeniusDiscriminant, FrobeniusTrace};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, FrobeniusTraceCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Frobenius-side case packaged from one finite short-Weierstrass curve.
#[derive(Clone)]
pub struct FrobeniusCurveCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub trace: FrobeniusTrace,
    pub discriminant: FrobeniusDiscriminant,
}

impl<F: Field> fmt::Debug for FrobeniusCurveCase<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_struct("FrobeniusCurveCase").finish()
    }
}

/// Returns a finite-field curve together with its Frobenius trace package.
pub fn arb_frobenius_curve_case<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<FrobeniusCurveCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>: CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
        + FrobeniusTraceCurveModel
        + 'static,
{
    arb_nonsingular_curve::<F>(config)
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
    let _ = |case: FrobeniusCurveCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.trace;
        let _ = case.discriminant;
    };
}
