use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, CurveModel, EnumerableCurveModel, ShortWeierstrassCurve,
};
use crate::fields::Fp;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// Returns a coupled short-Weierstrass curve together with one rational point.
pub fn arb_curve_and_point<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(ShortWeierstrassCurve<Fp<P>>, AffinePoint<Fp<P>>)> {
    arb_nonsingular_curve::<P>(config)
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
