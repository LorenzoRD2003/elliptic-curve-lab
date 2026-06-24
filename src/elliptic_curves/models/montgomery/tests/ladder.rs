use super::shared::{F5, F7, f5_curve, f7_scaled_curve, normalize_point};
use crate::elliptic_curves::{
    AffinePoint, MontgomeryXzPoint,
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use crate::fields::traits::Field;

fn xz_of<F: Field>(point: &AffinePoint<F>) -> MontgomeryXzPoint<F>
where
    F::Elem: Clone,
{
    MontgomeryXzPoint::from_affine_point(point)
}

#[test]
fn normalized_ladder_x_matches_affine_scalar_multiplication_on_a_small_example() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the normalized curve");

    for scalar in 0..8 {
        let expected = ambient
            .mul_scalar(&point, scalar)
            .expect("affine scalar multiplication should succeed");
        let actual = normalized.ladder_x(F5::from_i64(2), scalar);

        assert_eq!(actual, xz_of(&expected), "scalar = {scalar}");
    }
}

#[test]
fn normalized_ladder_pair_tracks_neighboring_multiples() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the normalized curve");

    for scalar in 0..8 {
        let (multiple, next_multiple) = normalized.ladder_xz_pair(F5::from_i64(2), scalar);
        let expected_multiple = ambient
            .mul_scalar(&point, scalar)
            .expect("affine scalar multiplication should succeed");
        let expected_next_multiple = ambient
            .mul_scalar(&point, scalar + 1)
            .expect("affine scalar multiplication should succeed");

        assert_eq!(multiple, xz_of(&expected_multiple), "scalar = {scalar}");
        assert_eq!(
            next_multiple,
            xz_of(&expected_next_multiple),
            "scalar = {scalar}"
        );
    }
}

#[test]
fn normalized_ladder_depends_only_on_the_x_coordinate_not_the_sign_of_y() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the normalized curve");
    let negated = ambient.neg(&point);

    for scalar in 0..8 {
        let from_point = normalized.ladder_x(F5::from_i64(2), scalar);
        let from_negated = match negated {
            AffinePoint::Finite { x, .. } => normalized.ladder_x(x, scalar),
            AffinePoint::Infinity => panic!("negation of a finite point should stay finite"),
        };

        assert_eq!(from_point, from_negated, "scalar = {scalar}");
    }
}

#[test]
fn normalized_ladder_returns_infinity_for_a_scalar_multiple_of_the_point_order() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the normalized curve");

    assert_eq!(
        ambient
            .mul_scalar(&point, 8)
            .expect("affine scalar multiplication should succeed"),
        ambient.identity()
    );

    assert_eq!(
        normalized.ladder_x(F5::from_i64(2), 8),
        MontgomeryXzPoint::Infinity
    );
}

#[test]
fn source_curve_try_ladder_x_matches_the_normalized_ladder_route() {
    let source = f7_scaled_curve();
    let normalized = source
        .try_as_normalized_montgomery()
        .expect("B = 2 is a square in F7");
    let source_point = source
        .point(F7::from_i64(2), F7::from_i64(2))
        .expect("sample point should lie on the source Montgomery curve");
    let normalized_point = normalize_point(&source, &normalized, &source_point)
        .expect("point should transport to the normalized target");

    let source_result = source
        .try_ladder_x(F7::from_i64(2), 3)
        .expect("ladder should be available when the same-field normalization exists");
    let normalized_result = normalized.ladder_x(F7::from_i64(2), 3);
    let ambient = normalized.as_montgomery_curve();
    let expected = ambient
        .mul_scalar(&normalized_point, 3)
        .expect("affine scalar multiplication should succeed");

    assert_eq!(source_result, normalized_result);
    assert_eq!(source_result, xz_of(&expected));
}
