use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;

use crate::elliptic_curves::traits::{AffineCurveModel, CurveModel, GroupCurveModel};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};
use crate::fields::{
    Fp,
    traits::{EnumerableFiniteField, Field, SqrtField},
};

pub(super) type F2 = Fp<2>;
pub(super) type F3 = Fp<3>;
pub(super) type F5 = Fp<5>;
pub(super) type F7 = Fp<7>;
pub(super) type F13 = Fp<13>;
pub(super) type F17 = Fp<17>;
pub(super) type F19 = Fp<19>;
pub(super) type F37 = Fp<37>;
pub(super) type F43 = Fp<43>;

pub(super) fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

pub(super) fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

pub(super) fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

pub(super) fn f43_curve() -> ShortWeierstrassCurve<F43> {
    ShortWeierstrassCurve::<F43>::new(F43::from_i64(2), F43::from_i64(3)).expect("valid curve")
}

pub(super) fn f13_j_1728_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::zero()).expect("valid curve")
}

pub(super) fn f13_j_zero_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(2)).expect("valid curve")
}

pub(super) fn f13_generic_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::from_i64(3)).expect("valid curve")
}

pub(super) fn f19_curve() -> ShortWeierstrassCurve<F19> {
    ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
}

pub(super) fn f37_curve() -> ShortWeierstrassCurve<F37> {
    ShortWeierstrassCurve::<F37>::new(F37::from_i64(2), F37::from_i64(3)).expect("valid curve")
}

pub(super) fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("point should lie on the curve")
}

pub(super) fn assert_contains(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    assert!(curve.contains(point));
}

pub(super) fn assert_group_law(
    curve: &ShortWeierstrassCurve<F7>,
    left: &AffinePoint<F7>,
    right: &AffinePoint<F7>,
    expected: &AffinePoint<F7>,
) {
    assert_contains(curve, left);
    assert_contains(curve, right);
    assert_contains(curve, expected);
    assert_eq!(curve.add(left, right), Ok(expected.clone()));
    assert_eq!(curve.sub(expected, right), Ok(left.clone()));
    assert_eq!(curve.sub(expected, left), Ok(right.clone()));
}

pub(super) fn assert_add_commutative(
    curve: &ShortWeierstrassCurve<F7>,
    left: &AffinePoint<F7>,
    right: &AffinePoint<F7>,
) {
    let left_right = curve
        .add(left, right)
        .expect("enumerated points should add successfully");
    let right_left = curve
        .add(right, left)
        .expect("enumerated points should add successfully");

    assert_eq!(left_right, right_left);
    assert_contains(curve, &left_right);
}

pub(super) fn assert_add_associative(
    curve: &ShortWeierstrassCurve<F7>,
    left: &AffinePoint<F7>,
    middle: &AffinePoint<F7>,
    right: &AffinePoint<F7>,
) {
    let left_grouped = curve
        .add(
            &curve
                .add(left, middle)
                .expect("enumerated points should add successfully"),
            right,
        )
        .expect("enumerated points should add successfully");
    let right_grouped = curve
        .add(
            left,
            &curve
                .add(middle, right)
                .expect("enumerated points should add successfully"),
        )
        .expect("enumerated points should add successfully");

    assert_eq!(left_grouped, right_grouped);
    assert_contains(curve, &left_grouped);
}

pub(super) fn assert_identity_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    assert_eq!(curve.add(&curve.identity(), point), Ok(point.clone()));
    assert_eq!(curve.add(point, &curve.identity()), Ok(point.clone()));
}

pub(super) fn assert_inverse_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    let inverse = curve.neg(point);

    assert_eq!(curve.add(point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.add(&inverse, point), Ok(curve.identity()));
}

pub(super) fn assert_scalar_mul_consistent(
    curve: &ShortWeierstrassCurve<F7>,
    point: &AffinePoint<F7>,
    n: u64,
    m: u64,
) {
    let left = curve
        .mul_scalar(point, n + m)
        .expect("scalar multiplication should succeed");
    let right = curve
        .add(
            &curve
                .mul_scalar(point, n)
                .expect("scalar multiplication should succeed"),
            &curve
                .mul_scalar(point, m)
                .expect("scalar multiplication should succeed"),
        )
        .expect("point addition should succeed");

    assert_eq!(left, right);
    assert_contains(curve, &left);
}

pub(super) fn first_nonsquare<F>() -> F::Elem
where
    F: EnumerableFiniteField + SqrtField,
{
    F::elements()
        .into_iter()
        .find(|value| !F::is_zero(value) && !F::has_square_root(value))
        .expect("small odd prime fields should contain non-squares")
}

pub(super) fn small_cyclic_distribution() -> BTreeMap<usize, usize> {
    BTreeMap::from([(1, 1), (2, 1), (3, 2), (6, 2)])
}
