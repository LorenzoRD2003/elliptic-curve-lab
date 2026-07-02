use crate::fields::traits::*;
use std::collections::HashMap;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::{ShortWeierstrassQuadraticTwist, TwistKind},
    traits::FiniteGroupCurveModel,
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

pub(super) type F43 = crate::fields::Fp43;
pub(super) type F241 = crate::fields::Fp241;

pub(super) fn f241_curve() -> ShortWeierstrassCurve<F241> {
    ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
        .expect("valid F241 curve")
}

pub(super) fn max_order_point_index<F>(curve: &ShortWeierstrassCurve<F>) -> usize
where
    F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField,
{
    curve
        .point_orders()
        .into_iter()
        .enumerate()
        .max_by_key(|(_, (_, order))| *order)
        .map(|(index, _)| index)
        .expect("small finite curve should have at least one point")
}

pub(super) fn genuine_twist_curve(
    curve: &ShortWeierstrassCurve<F241>,
) -> ShortWeierstrassCurve<F241> {
    F241::elements()
        .into_iter()
        .find_map(|candidate| {
            if F241::is_zero(&candidate) {
                return None;
            }
            let Ok(package) = ShortWeierstrassQuadraticTwist::new(curve.clone(), candidate) else {
                return None;
            };
            (package.kind() == TwistKind::Quadratic).then(|| package.twist().clone())
        })
        .expect("a prime-field curve should admit a genuine quadratic twist")
}

pub(super) fn sampler_covering_each_curve_by_index() -> impl FnMut(usize) -> Option<usize> {
    let mut next_index_by_upper_bound = HashMap::<usize, usize>::new();
    move |upper_bound: usize| {
        let next_index = next_index_by_upper_bound.entry(upper_bound).or_insert(0);
        let sampled = *next_index % upper_bound;
        *next_index += 1;
        Some(sampled)
    }
}
