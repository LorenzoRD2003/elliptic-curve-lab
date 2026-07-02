use super::shared::{F5, F7, f3_curve, f5_curve, f7_scaled_curve, normalize_point};
use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    montgomery::{MontgomeryXzPoint, NormalizedMontgomeryCurve},
    traits::{AffineCurveModel, CurveModel, EnumerableCurveModel, GroupCurveModel},
};
use crate::fields::traits::*;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_montgomery_curve_and_point,
};
use proptest::prelude::*;

fn xz_of<F: Field>(point: &AffinePoint<F>) -> MontgomeryXzPoint<F>
where
    F::Elem: Clone,
{
    MontgomeryXzPoint::from_affine_point(point)
}

fn normalized_case_f5() -> impl Strategy<
    Value = (
        MontgomeryCurve<F5>,
        NormalizedMontgomeryCurve<F5>,
        AffinePoint<F5>,
        AffinePoint<F5>,
        u64,
    ),
> {
    arb_montgomery_curve_and_point::<crate::fields::Fp5>(CurveStrategyConfig {
        include_identity_points: false,
        ..CurveStrategyConfig::default()
    })
    .prop_filter_map(
        "curve should admit same-field normalization",
        |(curve, point)| {
            let normalized = curve.try_as_normalized_montgomery().ok()?;
            let normalized_point = normalize_point(&curve, &normalized, &point).ok()?;
            Some((curve, normalized, point, normalized_point))
        },
    )
    .prop_flat_map(|(curve, normalized, point, normalized_point)| {
        (
            Just(curve),
            Just(normalized),
            Just(point),
            Just(normalized_point),
            0u64..16,
        )
    })
}

fn assert_exhaustive_ladder_agrees_with_affine_for_curve<F>(curve: MontgomeryCurve<F>)
where
    F: Field + SqrtField + EnumerableFiniteField,
    F::Elem: Clone,
{
    let normalized = curve
        .try_as_normalized_montgomery()
        .expect("test helper requires a same-field normalization");
    let ambient = normalized.as_montgomery_curve();

    for point in ambient.points() {
        if ambient.is_identity(&point) {
            continue;
        }

        let base_x = match &point {
            AffinePoint::Finite { x, .. } => x.clone(),
            AffinePoint::Infinity => continue,
        };
        let negated = ambient.neg(&point);
        let negated_x = match &negated {
            AffinePoint::Finite { x, .. } => x.clone(),
            AffinePoint::Infinity => continue,
        };

        for scalar in 0..=ambient.order() as u64 {
            let expected = ambient
                .mul_scalar(&point, scalar)
                .expect("affine scalar multiplication should succeed");
            let actual = normalized.ladder_x(base_x.clone(), scalar);
            let from_negated = normalized.ladder_x(negated_x.clone(), scalar);

            assert_eq!(
                actual,
                xz_of(&expected),
                "scalar = {scalar}, point = {point:?}"
            );
            assert_eq!(
                actual, from_negated,
                "scalar = {scalar}, point = {point:?}, negated = {negated:?}"
            );
        }
    }
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
fn exhaustive_ladder_checks_hold_over_f3() {
    assert_exhaustive_ladder_agrees_with_affine_for_curve(f3_curve());
}

#[test]
fn exhaustive_ladder_checks_hold_over_f5() {
    assert_exhaustive_ladder_agrees_with_affine_for_curve(f5_curve());
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
fn ladder_report_records_the_neighboring_pair_and_stays_honest_about_scope() {
    let normalized = f5_curve()
        .try_as_normalized_montgomery()
        .expect("B = 1 should normalize over the same field");
    let report = normalized.ladder_x_report(F5::from_i64(2), 3);

    assert!(F5::eq(report.base_x(), &F5::from_i64(2)));
    assert_eq!(report.scalar(), 3);
    assert_eq!(
        report.multiple_x(),
        &normalized.ladder_x(F5::from_i64(2), 3)
    );
    assert_eq!(
        report.next_multiple_x(),
        &normalized.ladder_xz_pair(F5::from_i64(2), 3).1
    );
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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_ladder_matches_affine_scalar_multiplication_on_random_normalizable_f5_curves(
        (curve, normalized, point, normalized_point, scalar) in normalized_case_f5(),
    ) {
        let base_x = match &point {
            AffinePoint::Finite { x, .. } => x.clone(),
            AffinePoint::Infinity => unreachable!("strategy excludes identity points"),
        };
        let ambient = normalized.as_montgomery_curve();
        let expected = ambient
            .mul_scalar(&normalized_point, scalar)
            .expect("affine scalar multiplication should succeed");

        prop_assert_eq!(curve.try_ladder_x(base_x, scalar).expect("normalizable curve should expose ladder"), xz_of(&expected));
    }

    #[test]
    fn property_ladder_is_sign_agnostic_on_random_normalizable_f5_curves(
        (_curve, normalized, _point, normalized_point, scalar) in normalized_case_f5(),
    ) {
        let ambient = normalized.as_montgomery_curve();
        let negated = ambient.neg(&normalized_point);
        let base_x = match &normalized_point {
            AffinePoint::Finite { x, .. } => x.clone(),
            AffinePoint::Infinity => unreachable!("strategy excludes identity points"),
        };
        let negated_x = match &negated {
            AffinePoint::Finite { x, .. } => x.clone(),
            AffinePoint::Infinity => unreachable!("negation of a finite point should stay finite"),
        };

        prop_assert_eq!(
            normalized.ladder_x(base_x, scalar),
            normalized.ladder_x(negated_x, scalar),
        );
    }
}
