use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, CurveError, ShortWeierstrassCurve};
use crate::fields::traits::Field;

use super::shared::{
    F7, F17, assert_add_associative, assert_add_commutative, assert_group_law, assert_identity_law,
    assert_inverse_law, assert_scalar_mul_consistent, f7_curve, f7_point, f43_curve,
};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};
use proptest::prelude::*;

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
