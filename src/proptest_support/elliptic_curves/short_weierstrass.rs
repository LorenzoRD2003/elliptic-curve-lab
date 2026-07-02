use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::EnumerableFiniteField;
use crate::proptest_support::config::CurveStrategyConfig;

/// Returns a non-singular short-Weierstrass curve over an enumerable field.
pub fn arb_nonsingular_curve<F>(
    _config: CurveStrategyConfig,
) -> BoxedStrategy<ShortWeierstrassCurve<F>>
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
            ShortWeierstrassCurve::<F>::new(a, b).ok()
        })
        .boxed()
}
