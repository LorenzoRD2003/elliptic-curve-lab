use super::shared::{F5, F7, f5_curve, f7_scaled_curve, normalize_point};
use crate::elliptic_curves::{
    AffinePoint,
    montgomery::{MontgomeryDifferentialArithmeticError, MontgomeryXzPoint},
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
fn x_dbl_matches_affine_doubling_on_the_normalized_curve() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the normalized curve");

    let doubled_x = normalized.x_dbl(&xz_of(&point));
    let expected = ambient
        .double(&point)
        .expect("affine doubling should succeed");

    assert_eq!(doubled_x, xz_of(&expected));
}

#[test]
fn x_dbl_sends_a_two_torsion_x_coordinate_to_infinity() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let two_torsion = ambient
        .point(F5::zero(), F5::zero())
        .expect("sample two-torsion point should lie on the curve");

    assert_eq!(
        normalized.x_dbl(&xz_of(&two_torsion)),
        xz_of(&ambient.identity())
    );
}

#[test]
fn x_add_matches_affine_addition_when_the_difference_x_is_known() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let left = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("left point should lie on the curve");
    let right = ambient
        .point(F5::from_i64(4), F5::from_i64(3))
        .expect("right point should lie on the curve");
    let difference = ambient
        .sub(&left, &right)
        .expect("difference point should be defined");
    let expected_sum = ambient.add(&left, &right).expect("sum should be defined");

    let actual = normalized
        .x_add(&xz_of(&left), &xz_of(&right), &xz_of(&difference))
        .expect("xADD should succeed when x(P-Q) is known");

    assert_eq!(actual, xz_of(&expected_sum));
}

#[test]
fn x_dbl_add_matches_separate_x_dbl_and_x_add_results() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let left = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("left point should lie on the curve");
    let right = ambient
        .point(F5::from_i64(4), F5::from_i64(3))
        .expect("right point should lie on the curve");
    let difference = ambient
        .sub(&left, &right)
        .expect("difference point should be defined");

    let (doubled, added) = normalized
        .x_dbl_add(&xz_of(&left), &xz_of(&right), &xz_of(&difference))
        .expect("xDBLADD should succeed when x(P-Q) is known");

    assert_eq!(doubled, normalized.x_dbl(&xz_of(&left)));
    assert_eq!(
        added,
        normalized
            .x_add(&xz_of(&left), &xz_of(&right), &xz_of(&difference))
            .expect("xADD should succeed when x(P-Q) is known")
    );
}

#[test]
fn x_add_handles_the_identity_special_case_when_the_difference_matches() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the curve");

    let sum = normalized
        .x_add(&xz_of(&ambient.identity()), &xz_of(&point), &xz_of(&point))
        .expect("O + P should be supported when x(P-O) = x(P) is supplied");

    assert_eq!(sum, xz_of(&point));
}

#[test]
fn x_add_rejects_incompatible_identity_difference_data() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let ambient = normalized.as_montgomery_curve();
    let point = ambient
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the curve");
    let wrong_difference = ambient
        .point(F5::zero(), F5::zero())
        .expect("sample point should lie on the curve");

    assert_eq!(
        normalized.x_add(
            &xz_of(&ambient.identity()),
            &xz_of(&point),
            &xz_of(&wrong_difference)
        ),
        Err(MontgomeryDifferentialArithmeticError::IncompatibleDifference)
    );
}

#[test]
fn differential_arithmetic_works_on_a_normalized_target_obtained_from_b_not_equal_one() {
    let source = f7_scaled_curve();
    let normalized = source
        .try_as_normalized_montgomery()
        .expect("B = 2 is a square in F7");
    let ambient = normalized.as_montgomery_curve();
    let source_point = source
        .point(F7::from_i64(2), F7::from_i64(2))
        .expect("sample point should lie on the source Montgomery curve");
    let normalized_point = normalize_point(&source, &normalized, &source_point)
        .expect("point should transport to the normalized target");
    let doubled_x = normalized.x_dbl(&xz_of(&normalized_point));
    let expected = ambient
        .double(&normalized_point)
        .expect("affine doubling should succeed on the normalized target");

    assert_eq!(doubled_x, xz_of(&expected));
}
