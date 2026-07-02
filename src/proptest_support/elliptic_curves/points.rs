use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Returns a coupled short-Weierstrass curve together with one rational point.
pub fn arb_curve_and_point<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(ShortWeierstrassCurve<F>, AffinePoint<F>)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(config)
        .prop_flat_map(move |curve| {
            let all_points = curve.points();
            let preferred_points = if config.include_identity_points {
                all_points.clone()
            } else {
                let finite_points = all_points
                    .iter()
                    .filter(|point| !curve.is_identity(point))
                    .cloned()
                    .collect::<Vec<_>>();
                if finite_points.is_empty() {
                    all_points.clone()
                } else {
                    finite_points
                }
            };

            (Just(curve), prop::sample::select(preferred_points))
        })
        .boxed()
}
