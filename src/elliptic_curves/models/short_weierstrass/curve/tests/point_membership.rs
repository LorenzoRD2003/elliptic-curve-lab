use crate::elliptic_curves::traits::{AffineCurveModel, CurveModel};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};
use crate::fields::traits::Field;

use super::shared::{F7, f7_curve};

#[test]
fn contains_accepts_affine_and_infinite_points_on_the_curve() {
    let curve = f7_curve();
    let affine_point = super::shared::f7_point(2, 1);
    let infinity = AffinePoint::<F7>::infinity();

    super::shared::assert_contains(&curve, &affine_point);
    super::shared::assert_contains(&curve, &infinity);
}

#[test]
fn contains_rejects_points_off_the_curve() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert!(!curve.contains(&point));
    assert!(!curve.is_on_curve_nonzero(&point));
}

#[test]
fn is_on_curve_nonzero_distinguishes_identity_from_finite_points() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let finite_point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let identity = AffinePoint::<F7>::infinity();

    assert!(curve.contains(&identity));
    assert!(curve.is_identity(&identity));
    assert!(!curve.is_on_curve_nonzero(&identity));

    assert!(curve.contains(&finite_point));
    assert!(!curve.is_identity(&finite_point));
    assert!(curve.is_on_curve_nonzero(&finite_point));
}
