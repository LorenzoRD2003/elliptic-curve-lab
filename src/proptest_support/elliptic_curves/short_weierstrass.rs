use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::{Field, Fp};
use crate::proptest_support::config::CurveStrategyConfig;

/// Returns a non-singular short-Weierstrass curve over `GF(P)`.
pub fn arb_nonsingular_curve<const P: u64>(
    _config: CurveStrategyConfig,
) -> BoxedStrategy<ShortWeierstrassCurve<Fp<P>>> {
    (0..P, 0..P)
        .prop_filter_map("curve must be non-singular", |(a, b)| {
            ShortWeierstrassCurve::<Fp<P>>::new(
                Fp::<P>::elem_from_u64(a),
                Fp::<P>::elem_from_u64(b),
            )
            .ok()
        })
        .boxed()
}
