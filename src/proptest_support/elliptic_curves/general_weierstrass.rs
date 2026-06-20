use proptest::prelude::*;

use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};
use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve};
use crate::fields::{Fp, traits::Field};
use crate::proptest_support::config::CurveStrategyConfig;

/// Returns a non-singular general Weierstrass curve over `GF(P)`.
pub fn arb_nonsingular_general_weierstrass_curve<const P: u64>(
    _config: CurveStrategyConfig,
) -> BoxedStrategy<GeneralWeierstrassCurve<Fp<P>>> {
    (0..P, 0..P, 0..P, 0..P, 0..P)
        .prop_filter_map("curve must be non-singular", |(a1, a2, a3, a4, a6)| {
            GeneralWeierstrassCurve::<Fp<P>>::new(
                Fp::<P>::elem_from_u64(a1),
                Fp::<P>::elem_from_u64(a2),
                Fp::<P>::elem_from_u64(a3),
                Fp::<P>::elem_from_u64(a4),
                Fp::<P>::elem_from_u64(a6),
            )
            .ok()
        })
        .boxed()
}

/// Returns a coupled general-Weierstrass curve together with one rational
/// point obtained by exhaustive lifting over `GF(P)`.
pub fn arb_general_weierstrass_curve_and_point<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(GeneralWeierstrassCurve<Fp<P>>, AffinePoint<Fp<P>>)> {
    arb_nonsingular_general_weierstrass_curve::<P>(config)
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
