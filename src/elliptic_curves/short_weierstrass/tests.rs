use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveError, CurveIsomorphism, CurveIsomorphismError, CurveModel,
    EnumerableCurveModel, FiniteAbelianGroupStructure, FiniteGroupCurveModel,
    FrobeniusTraceCurveModel, GroupCurveModel, HasJInvariant, LiftXCoordinate, PointCountStrategy,
};
use crate::fields::{EnumerableFiniteField, Field, Fp, Q, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::arb_nonsingular_curve;

type F2 = Fp<2>;
type F3 = Fp<3>;
type F5 = Fp<5>;
type F7 = Fp<7>;
type F13 = Fp<13>;
type F17 = Fp<17>;
type F19 = Fp<19>;
type F37 = Fp<37>;
type F43 = Fp<43>;

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

fn f5_noncyclic_curve() -> ShortWeierstrassCurve<F5> {
    ShortWeierstrassCurve::<F5>::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
}

fn f43_curve() -> ShortWeierstrassCurve<F43> {
    ShortWeierstrassCurve::<F43>::new(F43::from_i64(2), F43::from_i64(3)).expect("valid curve")
}

fn f13_j_1728_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::zero()).expect("valid curve")
}

fn f13_j_zero_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(2)).expect("valid curve")
}

fn f13_generic_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::from_i64(2), F13::from_i64(3)).expect("valid curve")
}

fn f19_curve() -> ShortWeierstrassCurve<F19> {
    ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
}

fn f37_curve() -> ShortWeierstrassCurve<F37> {
    ShortWeierstrassCurve::<F37>::new(F37::from_i64(2), F37::from_i64(3)).expect("valid curve")
}

fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("point should lie on the curve")
}

fn assert_contains(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    assert!(curve.contains(point));
}

fn assert_group_law(
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

fn assert_add_commutative(
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

fn assert_add_associative(
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

fn assert_identity_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    assert_eq!(curve.add(&curve.identity(), point), Ok(point.clone()));
    assert_eq!(curve.add(point, &curve.identity()), Ok(point.clone()));
}

fn assert_inverse_law(curve: &ShortWeierstrassCurve<F7>, point: &AffinePoint<F7>) {
    let inverse = curve.neg(point);

    assert_eq!(curve.add(point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.add(&inverse, point), Ok(curve.identity()));
}

fn assert_scalar_mul_consistent(
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

#[test]
fn constructor_rejects_characteristics_two_and_three() {
    assert!(matches!(
        ShortWeierstrassCurve::<F2>::new(F2::zero(), F2::one()),
        Err(CurveError::UnsupportedCharacteristic { characteristic: 2 }),
    ));
    assert!(matches!(
        ShortWeierstrassCurve::<F3>::new(F3::zero(), F3::one()),
        Err(CurveError::UnsupportedCharacteristic { characteristic: 3 }),
    ));
}

#[test]
fn constructor_rejects_singular_coefficients() {
    assert!(matches!(
        ShortWeierstrassCurve::<F5>::new(F5::zero(), F5::zero()),
        Err(CurveError::SingularCurve),
    ));
}

#[test]
fn accessors_discriminant_and_rhs_match_the_model() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(F7::eq(curve.a(), &F7::from_i64(2)));
    assert!(F7::eq(curve.b(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c4(), &F7::from_i64(2)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(5)));
    assert!(F7::eq(
        &LiftXCoordinate::rhs(&curve, &F7::from_i64(2)),
        &F7::from_i64(1)
    ));
}

#[test]
fn weierstrass_invariants_satisfy_the_classical_relation() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
    let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

    assert!(F7::eq(&left, &right));
}

#[test]
fn j_invariant_matches_a_classical_exact_example_over_q() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

    assert!(Q::eq(&curve.discriminant(), &q(64, 1)));
    assert!(Q::eq(&curve.c4(), &q(48, 1)));
    assert!(Q::eq(&curve.c6(), &q(0, 1)));
    assert!(Q::eq(&curve.j_invariant(), &q(1728, 1)));
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_method() {
    let curve = f7_curve();

    assert!(F7::eq(
        &HasJInvariant::j_invariant(&curve),
        &ShortWeierstrassCurve::j_invariant(&curve),
    ));
}

#[test]
fn scaled_by_rejects_noninvertible_scale() {
    let curve = f7_curve();

    assert!(matches!(
        curve.scaled_by(F7::zero()),
        Err(CurveIsomorphismError::NonInvertibleScale)
    ));
}

#[test]
fn scaled_by_applies_the_expected_u4_and_u6_coefficients() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F7::eq(scaled.a(), &F7::from_i64(1)));
    assert!(F7::eq(scaled.b(), &F7::from_i64(3)));
}

#[test]
fn isomorphic_via_scale_matches_the_scaled_curve() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.isomorphic_via_scale(&scaled, &F7::from_i64(3)));
    assert!(!curve.isomorphic_via_scale(&scaled, &F7::from_i64(2)));
}

#[test]
fn isomorphic_via_scale_returns_false_for_noninvertible_scale() {
    let curve = f7_curve();

    assert!(!curve.isomorphic_via_scale(&curve, &F7::zero()));
}

fn first_nonsquare<F>() -> F::Elem
where
    F: EnumerableFiniteField + SqrtField,
{
    F::elements()
        .into_iter()
        .find(|value| !F::is_zero(value) && !F::has_square_root(value))
        .expect("small odd prime fields should contain non-squares")
}

#[test]
fn quadratic_twist_rejects_zero_factor() {
    let curve = f19_curve();

    assert!(matches!(
        curve.quadratic_twist(F19::zero()),
        Err(CurveIsomorphismError::NonInvertibleScale)
    ));
}

#[test]
fn quadratic_twist_preserves_the_j_invariant_over_f19() {
    let curve = f19_curve();
    let twist = curve
        .quadratic_twist(F19::from_i64(2))
        .expect("non-zero twist factor should produce a valid model");

    assert!(curve.has_same_j_invariant(&twist));
}

#[test]
fn quadratic_twist_by_a_square_is_base_field_isomorphic_over_f19() {
    let curve = f19_curve();
    let square = F19::from_i64(4);
    let twist = curve
        .quadratic_twist(square)
        .expect("square twist factor should produce a valid model");

    assert!(F19::has_square_root(&square));
    assert!(curve.is_isomorphic_to(&twist));
}

#[test]
fn quadratic_twist_by_a_nonsquare_is_not_base_field_isomorphic_in_the_sample_f19_case() {
    let curve = f19_curve();
    let nonsquare = first_nonsquare::<F19>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert!(!F19::has_square_root(&nonsquare));
    assert!(curve.has_same_j_invariant(&twist));
    assert!(!curve.is_isomorphic_to(&twist));
}

#[test]
fn quadratic_twist_point_count_relation_holds_over_f19() {
    let curve = f19_curve();
    let nonsquare = first_nonsquare::<F19>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert_eq!(curve.order() + twist.order(), 2 * 19 + 2);
}

#[test]
fn quadratic_twist_point_count_relation_holds_over_f37() {
    let curve = f37_curve();
    let nonsquare = first_nonsquare::<F37>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert_eq!(curve.order() + twist.order(), 2 * 37 + 2);
}

#[test]
fn find_isomorphism_to_recovers_a_base_field_scaling_witness() {
    let curve = f7_curve();
    let other = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");
    let point = f7_point(2, 1);

    let isomorphism = curve
        .find_isomorphism_to(&other)
        .expect("a base-field scaling witness should exist");
    let image = isomorphism
        .evaluate(&point)
        .expect("the witness isomorphism should transport domain points");

    assert!(isomorphism.codomain().contains(&image));
    assert_eq!(
        isomorphism
            .inverse()
            .expect("inverse should exist")
            .evaluate(&image)
            .expect("inverse should recover the original point"),
        point
    );
}

#[test]
fn find_isomorphism_to_returns_none_for_same_j_but_no_base_field_scale() {
    let curve = f7_curve();
    let same_j_not_base_isomorphic =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4)).expect("valid curve");

    assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
    assert!(
        curve
            .find_isomorphism_to(&same_j_not_base_isomorphic)
            .is_none()
    );
}

#[test]
fn is_isomorphic_to_returns_true_when_a_base_field_witness_exists() {
    let curve = f7_curve();
    let other = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.is_isomorphic_to(&other));
}

#[test]
fn is_isomorphic_to_returns_false_when_only_the_j_invariant_matches() {
    let curve = f7_curve();
    let same_j_not_base_isomorphic =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4)).expect("valid curve");

    assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
    assert!(!curve.is_isomorphic_to(&same_j_not_base_isomorphic));
}

#[test]
fn generic_curve_has_only_plus_minus_one_automorphisms_over_f7() {
    let curve = f7_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 2);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F7::eq(iso.scaling_factor(), &F7::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F7::eq(iso.scaling_factor(), &F7::from_i64(-1)))
    );
}

#[test]
fn special_j_1728_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_j_1728_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F13::is_zero(curve.b()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_1728_curve_has_four_automorphisms_over_f13() {
    let curve = f13_j_1728_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 4);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(5)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-5)))
    );
}

#[test]
fn special_j_zero_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_j_zero_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F13::is_zero(curve.a()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_zero_curve_has_six_automorphisms_over_f13() {
    let curve = f13_j_zero_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 6);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(4)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-4)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(3)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-3)))
    );
}

#[test]
fn generic_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_generic_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(!F13::is_zero(curve.a()));
    assert!(!F13::is_zero(curve.b()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_1728_same_j_does_not_force_base_field_isomorphism() {
    let first = f13_j_1728_curve();
    let second =
        ShortWeierstrassCurve::<F13>::new(F13::from_i64(1), F13::zero()).expect("valid curve");

    assert!(first.has_same_j_invariant(&second));
    assert!(!first.is_isomorphic_to(&second));
    assert!(first.find_isomorphism_to(&second).is_none());
}

#[test]
fn special_j_zero_same_j_does_not_force_base_field_isomorphism() {
    let first = f13_j_zero_curve();
    let second =
        ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(1)).expect("valid curve");

    assert!(first.has_same_j_invariant(&second));
    assert!(!first.is_isomorphic_to(&second));
    assert!(first.find_isomorphism_to(&second).is_none());
}

#[test]
fn has_same_j_invariant_detects_scaled_models() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.has_same_j_invariant(&scaled));
}

#[test]
fn has_same_j_invariant_distinguishes_curves_with_different_j() {
    let first = f7_curve();
    let second =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::from_i64(1)).expect("valid curve");

    assert!(!first.has_same_j_invariant(&second));
}

#[test]
fn contains_accepts_affine_and_infinite_points_on_the_curve() {
    let curve = f7_curve();
    let affine_point = f7_point(2, 1);
    let infinity = AffinePoint::<F7>::infinity();

    assert_contains(&curve, &affine_point);
    assert_contains(&curve, &infinity);
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
fn point_constructor_accepts_valid_affine_coordinates() {
    let curve = f7_curve();

    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");

    assert!(matches!(point, AffinePoint::Finite { .. }));
}

#[test]
fn point_constructor_rejects_invalid_affine_coordinates() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(matches!(
        curve.point(F7::from_i64(2), F7::from_i64(2)),
        Err(CurveError::PointNotOnCurve)
    ));
}

#[test]
fn characteristic_zero_fields_are_allowed() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");
    let point = curve
        .point(q(0, 1), q(0, 1))
        .expect("point should lie on the curve");

    assert!(curve.contains(&point));
    assert_eq!(Q::characteristic(), 0);
}

#[test]
fn point_from_x_returns_one_point_when_rhs_has_a_square_root() {
    let curve = f7_curve();

    let point = curve
        .point_from_x(F7::from_i64(2))
        .expect("x = 2 should lift to a point");

    assert_contains(&curve, &point);
    assert!(matches!(point, AffinePoint::Finite { .. }));
}

#[test]
fn point_from_x_returns_none_when_rhs_is_not_a_square() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(curve.point_from_x(F7::from_i64(0)).is_none());
}

#[test]
fn points_from_x_returns_both_points_when_they_are_distinct() {
    let curve = f7_curve();

    let (left, right) = curve
        .points_from_x(F7::from_i64(2))
        .expect("x = 2 should lift to two points");

    assert_contains(&curve, &left);
    assert_contains(&curve, &right);
    assert_ne!(left, right);
}

#[test]
fn points_from_x_repeats_the_point_when_the_square_root_is_zero() {
    let curve = f7_curve();

    let (left, right) = curve
        .points_from_x(F7::from_i64(6))
        .expect("x = 6 should give y = 0");

    assert_eq!(left, right);
    assert_contains(&curve, &left);
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

#[test]
fn points_from_x_works_over_q_when_an_exact_root_exists() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

    let (left, right) = curve
        .points_from_x(q(1, 1))
        .expect("x = 1 should give y = 0 in Q");

    assert_eq!(left, right);
    assert!(curve.contains(&left));
}

#[test]
fn finite_point_enumeration_lists_exactly_the_non_identity_points() {
    let curve = f7_curve();
    let finite_points = curve.finite_points();

    assert_eq!(finite_points.len(), 5);
    assert!(
        finite_points
            .iter()
            .all(|point| curve.is_on_curve_nonzero(point))
    );
}

#[test]
fn full_point_enumeration_includes_identity_and_order() {
    let curve = f7_curve();
    let points = curve.points();

    assert_eq!(points.len(), 6);
    assert!(curve.is_identity(points.first().expect("identity should be present")));
    assert_eq!(curve.order(), 6);
}

#[test]
fn public_point_count_api_prefers_character_sum_in_auto_mode() {
    let curve = f43_curve();

    let report = curve
        .count_points(PointCountStrategy::Auto)
        .expect("automatic point count should succeed");

    assert_eq!(report.strategy(), PointCountStrategy::QuadraticCharacter);
    assert_eq!(report.curve_order(), curve.order() as u128);
}

#[test]
fn public_frobenius_trace_by_agrees_with_the_exhaustive_trace() {
    let curve = f43_curve();

    assert_eq!(
        curve.frobenius_trace_by(PointCountStrategy::Exhaustive),
        curve.frobenius_trace()
    );
}

#[test]
fn random_point_uses_the_supplied_index_sampler() {
    let curve = f7_curve();
    let expected = curve.points()[2].clone();
    let mut sampler = |upper_bound: usize| {
        assert_eq!(upper_bound, 6);
        Some(2)
    };

    let sampled = curve
        .random_point(&mut sampler)
        .expect("sampler should choose an existing point");

    assert_eq!(sampled, expected);
}

#[test]
fn random_point_propagates_sampler_failure() {
    let curve = f7_curve();
    let mut sampler = |_upper_bound: usize| None;

    assert!(curve.random_point(&mut sampler).is_none());
}

#[test]
fn group_negation_matches_affine_involution() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert_eq!(curve.neg(&point), f7_point(2, 6));
    assert_eq!(curve.neg(&curve.identity()), curve.identity());
}

#[test]
fn group_add_handles_identity_and_inverse_cases() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert_identity_law(&curve, &point);
    assert_inverse_law(&curve, &point);
}

#[test]
fn group_add_and_double_match_known_small_field_examples() {
    let curve = f7_curve();
    let p = f7_point(2, 1);
    let q = f7_point(3, 1);
    let two_p = f7_point(3, 6);
    let p_plus_q = f7_point(2, 6);
    let torsion_point = f7_point(6, 0);

    assert_eq!(curve.double(&p), Ok(two_p));
    assert_group_law(&curve, &p, &q, &p_plus_q);
    assert_eq!(curve.sub(&p, &q), Ok(torsion_point));
}

#[test]
fn doubling_a_two_torsion_point_returns_the_identity() {
    let curve = f7_curve();
    let point = f7_point(6, 0);

    assert_eq!(curve.double(&point), Ok(curve.identity()));
}

#[test]
fn scalar_multiplication_matches_repeated_addition_examples() {
    let curve = f7_curve();
    let point = f7_point(2, 1);
    let three_p = f7_point(6, 0);
    let minus_two_p = f7_point(3, 1);

    assert_eq!(curve.mul_scalar(&point, 0), Ok(curve.identity()));
    assert_eq!(curve.mul_scalar(&point, 1), Ok(point.clone()));
    assert_eq!(curve.mul_scalar(&point, 3), Ok(three_p));
    assert_eq!(curve.mul_scalar(&point, 6), Ok(curve.identity()));
    assert_eq!(curve.mul_scalar_signed(&point, -2), Ok(minus_two_p));
    assert_scalar_mul_consistent(&curve, &point, 2, 3);
    assert_scalar_mul_consistent(&curve, &point, 1, 5);
}

#[test]
fn group_operations_reject_points_outside_the_curve() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let valid = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert_eq!(
        curve.add(&valid, &invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(curve.double(&invalid), Err(CurveError::PointNotOnCurve));
    assert_eq!(
        curve.sub(&valid, &invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        curve.mul_scalar(&invalid, 3),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        curve.mul_scalar_signed(&invalid, -3),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn enumerated_points_form_an_abelian_group_in_the_small_example() {
    let curve = f7_curve();
    let points = curve.points();

    for left in &points {
        for right in &points {
            assert_add_commutative(&curve, left, right);

            for third in &points {
                assert_add_associative(&curve, left, right, third);
            }
        }
    }
}

#[test]
fn torsion_helper_detects_known_orders_in_the_small_example() {
    let curve = f7_curve();
    let order_six_point = f7_point(2, 1);
    let order_two_point = f7_point(6, 0);
    let identity = curve.identity();

    assert!(curve.is_torsion_point(&order_six_point, 6));
    assert!(!curve.is_torsion_point(&order_six_point, 3));
    assert!(curve.is_torsion_point(&order_two_point, 2));
    assert!(curve.is_torsion_point(&identity, 5));
}

#[test]
fn torsion_helper_rejects_zero_order_and_invalid_points() {
    let curve = f7_curve();
    let valid = f7_point(2, 1);
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert!(!curve.is_torsion_point(&valid, 0));
    assert!(!curve.is_torsion_point(&invalid, 6));
}

#[test]
fn point_order_matches_known_small_group_examples() {
    let curve = f7_curve();
    let order_six_point = f7_point(2, 1);
    let order_two_point = f7_point(6, 0);

    assert_eq!(curve.point_order(&curve.identity()), Some(1));
    assert_eq!(curve.point_order(&order_two_point), Some(2));
    assert_eq!(curve.point_order(&order_six_point), Some(6));
}

#[test]
fn point_order_returns_none_for_points_outside_the_curve() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert_eq!(curve.point_order(&invalid), None);
}

#[test]
fn point_orders_cover_the_full_small_curve_group() {
    let curve = f7_curve();
    let point_orders = curve.point_orders();

    assert_eq!(point_orders.len(), curve.order());
    assert_eq!(point_orders[0], (curve.identity(), 1));
    assert!(point_orders.contains(&(f7_point(6, 0), 2)));
    assert!(point_orders.contains(&(f7_point(2, 1), 6)));
    assert!(point_orders.contains(&(f7_point(2, 6), 6)));
}

#[test]
fn points_of_order_filters_exact_orders() {
    let curve = f7_curve();

    assert_eq!(curve.points_of_order(1), vec![curve.identity()]);
    assert_eq!(curve.points_of_order(2), vec![f7_point(6, 0)]);
    assert_eq!(
        curve.points_of_order(3),
        vec![f7_point(3, 1), f7_point(3, 6)]
    );
    assert_eq!(
        curve.points_of_order(6),
        vec![f7_point(2, 1), f7_point(2, 6)]
    );
    assert!(curve.points_of_order(4).is_empty());
}

#[test]
fn order_distribution_matches_the_small_cyclic_example() {
    let curve = f7_curve();
    let expected = BTreeMap::from([(1, 1), (2, 1), (3, 2), (6, 2)]);

    assert_eq!(curve.order_distribution(), expected);
}

#[test]
fn exponent_generator_and_cyclicity_match_the_small_example() {
    let curve = f7_curve();
    let generator = curve.generator().expect("group should be cyclic");
    let structure = curve.group_structure();

    assert_eq!(curve.exponent(), 6);
    assert!(curve.is_cyclic());
    assert_eq!(
        structure,
        FiniteAbelianGroupStructure {
            order: 6,
            exponent: 6,
            cyclic: true,
            invariant_factors: None,
        }
    );
    assert_eq!(curve.describe_group_structure(), "Z/6Z");
    assert_eq!(curve.point_order(&generator), Some(curve.order()));
}

#[test]
fn noncyclic_f5_example_has_split_two_torsion_structure() {
    let curve = f5_noncyclic_curve();
    let expected = BTreeMap::from([(1, 1), (2, 3), (4, 4)]);
    let structure = curve.group_structure();

    assert_eq!(curve.order(), 8);
    assert_eq!(curve.order_distribution(), expected);
    assert_eq!(curve.exponent(), 4);
    assert_eq!(curve.generator(), None);
    assert!(!curve.is_cyclic());
    assert_eq!(
        structure,
        FiniteAbelianGroupStructure {
            order: 8,
            exponent: 4,
            cyclic: false,
            invariant_factors: Some((2, 4)),
        }
    );
    assert_eq!(curve.describe_group_structure(), "Z/2Z x Z/4Z");
}

#[test]
fn exhaustive_group_axiom_check_passes_for_a_small_f43_curve() {
    let curve = f43_curve();

    assert_eq!(curve.check_group_axioms(), Ok(()));
}

fn curve_and_group_data() -> impl Strategy<
    Value = (
        ShortWeierstrassCurve<F17>,
        AffinePoint<F17>,
        AffinePoint<F17>,
        u64,
        u64,
    ),
> {
    arb_nonsingular_curve::<17>(CurveStrategyConfig::default()).prop_flat_map(|curve| {
        let points = curve.points();
        let len = points.len();

        (
            Just(curve.clone()),
            Just(points),
            0usize..len,
            0usize..len,
            0u64..8,
            0u64..8,
        )
            .prop_map(|(curve, points, left_index, right_index, n, m)| {
                (
                    curve,
                    points[left_index].clone(),
                    points[right_index].clone(),
                    n,
                    m,
                )
            })
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_short_weierstrass_invariants_satisfy_the_classical_relation(
        curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default()),
    ) {
        let left = F17::sub(&F17::cube(&curve.c4()), &F17::square(&curve.c6()));
        let right = F17::mul(&F17::from_i64(1728), &curve.discriminant());

        prop_assert!(F17::eq(&left, &right));
    }

    #[test]
    fn property_short_weierstrass_group_law_holds_on_enumerated_points(
        (curve, left, right, n, m) in curve_and_group_data(),
    ) {
        let left_plus_right = curve.add(&left, &right).expect("enumerated points should add");
        let right_plus_left = curve.add(&right, &left).expect("enumerated points should add");
        let inverse = curve.neg(&left);
        let scalar_sum = curve.mul_scalar(&left, n + m).expect("scalar multiplication should succeed");
        let split_scalar = curve
            .add(
                &curve.mul_scalar(&left, n).expect("scalar multiplication should succeed"),
                &curve.mul_scalar(&left, m).expect("scalar multiplication should succeed"),
            )
            .expect("point addition should succeed");

        prop_assert_eq!(left_plus_right, right_plus_left);
        prop_assert_eq!(curve.add(&left, &inverse).expect("inverse sum should succeed"), curve.identity());
        prop_assert_eq!(scalar_sum, split_scalar);
    }
}
