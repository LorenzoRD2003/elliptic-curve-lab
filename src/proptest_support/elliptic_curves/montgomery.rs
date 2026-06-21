use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    traits::{CurveModel, EnumerableCurveModel},
};
use crate::fields::{Fp, traits::Field};
use crate::proptest_support::config::CurveStrategyConfig;

/// Returns a non-singular Montgomery curve over `GF(P)`.
pub fn arb_nonsingular_montgomery_curve<const P: u64>(
    _config: CurveStrategyConfig,
) -> BoxedStrategy<MontgomeryCurve<Fp<P>>> {
    (0..P, 0..P)
        .prop_filter_map("curve must be non-singular", |(a, b)| {
            MontgomeryCurve::<Fp<P>>::new(Fp::<P>::elem_from_u64(a), Fp::<P>::elem_from_u64(b)).ok()
        })
        .boxed()
}

/// Returns a coupled Montgomery curve together with one rational point
/// obtained by exhaustive lifting over `GF(P)`.
pub fn arb_montgomery_curve_and_point<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(MontgomeryCurve<Fp<P>>, AffinePoint<Fp<P>>)> {
    arb_nonsingular_montgomery_curve::<P>(config)
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
