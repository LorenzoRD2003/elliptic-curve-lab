use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    traits::{CurveModel, EnumerableCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;

/// Returns a non-singular Montgomery curve over an enumerable field.
pub fn arb_nonsingular_montgomery_curve<F>(
    _config: CurveStrategyConfig,
) -> BoxedStrategy<MontgomeryCurve<F>>
where
    F: EnumerableFiniteField + 'static,
    F::Elem: 'static,
{
    let elements = F::elements();
    (
        prop::sample::select(elements.clone()),
        prop::sample::select(elements),
    )
        .prop_filter_map("curve must be non-singular", |(a, b)| {
            MontgomeryCurve::<F>::new(a, b).ok()
        })
        .boxed()
}

/// Returns a coupled Montgomery curve together with one rational point
/// obtained by exhaustive lifting over an enumerable field.
pub fn arb_montgomery_curve_and_point<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(MontgomeryCurve<F>, AffinePoint<F>)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    MontgomeryCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_montgomery_curve::<F>(config)
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
