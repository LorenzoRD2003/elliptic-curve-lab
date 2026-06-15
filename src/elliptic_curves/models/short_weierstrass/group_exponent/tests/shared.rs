use num_bigint::BigUint;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    traits::{AffineCurveModel, EnumerableCurveModel, PointIndexSampler},
};
use crate::fields::{
    Fp,
    traits::{EnumerableFiniteField, Field, FiniteField, QuadraticCharacterFiniteField, SqrtField},
};

pub(super) type F7 = Fp<7>;

pub(super) fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

pub(super) fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

pub(super) fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("point should lie on the curve")
}

pub(super) fn alternate_f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::from_i64(1)).expect("valid curve")
}

pub(super) fn alternate_f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    alternate_f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("point should lie on the curve")
}

pub(super) fn sampler_from_indices(indices: Vec<usize>) -> impl PointIndexSampler {
    let mut indices = indices.into_iter();
    move |upper_bound: usize| indices.next().filter(|index| *index < upper_bound)
}

pub(super) fn point_index<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> usize
where
    F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField,
{
    curve
        .points()
        .iter()
        .position(|candidate| candidate == point)
        .expect("sample point should appear in the enumerated group")
}
