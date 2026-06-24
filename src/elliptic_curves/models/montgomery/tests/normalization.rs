use super::shared::{F5, F7, f5_curve, f7_nonsquare_scaled_curve, f7_scaled_curve};
use crate::elliptic_curves::{AffinePoint, CurveError, MontgomeryNormalizationError};
use crate::fields::traits::Field;

#[test]
fn normalization_of_b_equal_one_curve_keeps_the_same_a_and_produces_a_square_root_of_b() {
    let curve = f5_curve();
    let public_normalized = curve
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");

    let normalization = curve
        .try_normalize()
        .expect("B = 1 should normalize over the same field");

    assert_eq!(public_normalized, normalization.target().clone());
    assert_eq!(normalization.source(), &curve);
    assert!(F5::eq(normalization.target().a(), curve.a()));
    assert!(F5::eq(&F5::square(normalization.sqrt_b()), curve.b()));
}

#[test]
fn normalization_of_scaled_curve_records_a_same_field_square_root_of_b() {
    let curve = f7_scaled_curve();
    let public_normalized = curve
        .try_as_normalized_montgomery()
        .expect("B = 2 is a square in F7");

    let normalization = curve.try_normalize().expect("B = 2 is a square in F7");

    assert_eq!(public_normalized, normalization.target().clone());
    assert_eq!(normalization.source(), &curve);
    assert!(F7::eq(normalization.target().a(), curve.a()));
    assert!(F7::eq(&F7::square(normalization.sqrt_b()), curve.b()));
}

#[test]
fn normalization_transports_points_by_scaling_only_the_y_coordinate() {
    let curve = f5_curve();
    let normalization = curve
        .try_normalize()
        .expect("B = 1 should normalize over the same field");
    let source_point = AffinePoint::<F5>::new(F5::from_i64(2), F5::from_i64(2));

    let target_point = normalization
        .map_source_point(&source_point)
        .expect("source point should transport to the normalized target");

    match target_point {
        AffinePoint::Infinity => panic!("finite point should stay finite under normalization"),
        AffinePoint::Finite { x, y } => {
            assert!(F5::eq(&x, &F5::from_i64(2)));
            assert!(F5::eq(
                &y,
                &F5::mul(normalization.sqrt_b(), &F5::from_i64(2))
            ));
        }
    }

    let roundtrip = normalization
        .map_target_point(&target_point)
        .expect("target point should transport back to the source model");
    assert_eq!(roundtrip, source_point);
}

#[test]
fn normalization_rejects_non_square_b_honestly() {
    let curve = f7_nonsquare_scaled_curve();

    assert!(matches!(
        curve.try_as_normalized_montgomery(),
        Err(MontgomeryNormalizationError::NoSameFieldNormalization)
    ));
}

#[test]
fn normalization_reports_invalid_source_and_target_points_honestly() {
    let curve = f5_curve();
    let normalization = curve
        .try_normalize()
        .expect("B = 1 should normalize over the same field");
    let bad_source_point = AffinePoint::<F5>::new(F5::zero(), F5::one());
    let bad_target_point = AffinePoint::<F5>::new(F5::one(), F5::one());

    assert_eq!(
        normalization.map_source_point(&bad_source_point),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        normalization.map_target_point(&bad_target_point),
        Err(CurveError::PointNotOnCurve)
    );
}
