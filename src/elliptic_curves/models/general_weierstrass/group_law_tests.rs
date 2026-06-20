use proptest::prelude::*;

use crate::elliptic_curves::traits::{
    CurveModel, CurveModelConversion, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve};
use crate::fields::{Fp, traits::Field};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::arb_nonsingular_general_weierstrass_curve;

fn exhaustive_nonsingular_curves<const P: u64>() -> Vec<GeneralWeierstrassCurve<Fp<P>>> {
    let mut curves = Vec::new();

    for a1 in 0..P {
        for a2 in 0..P {
            for a3 in 0..P {
                for a4 in 0..P {
                    for a6 in 0..P {
                        if let Ok(curve) = GeneralWeierstrassCurve::<Fp<P>>::new(
                            Fp::<P>::elem_from_u64(a1),
                            Fp::<P>::elem_from_u64(a2),
                            Fp::<P>::elem_from_u64(a3),
                            Fp::<P>::elem_from_u64(a4),
                            Fp::<P>::elem_from_u64(a6),
                        ) {
                            curves.push(curve);
                        }
                    }
                }
            }
        }
    }

    curves
}

fn group_law_case<const P: u64>() -> impl Strategy<
    Value = (
        GeneralWeierstrassCurve<Fp<P>>,
        AffinePoint<Fp<P>>,
        AffinePoint<Fp<P>>,
        u64,
        u64,
    ),
> {
    arb_nonsingular_general_weierstrass_curve::<P>(CurveStrategyConfig::default()).prop_flat_map(
        |curve| {
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
        },
    )
}

#[test]
fn exhaustive_group_axiom_check_passes_for_every_nonsingular_curve_over_f2() {
    for curve in exhaustive_nonsingular_curves::<2>() {
        assert_eq!(curve.check_group_axioms(), Ok(()));
    }
}

#[test]
fn exhaustive_group_axiom_check_passes_for_every_nonsingular_curve_over_f3() {
    for curve in exhaustive_nonsingular_curves::<3>() {
        assert_eq!(curve.check_group_axioms(), Ok(()));
    }
}

#[test]
fn characteristic_two_point_orders_are_defined_for_every_enumerated_point() {
    for curve in exhaustive_nonsingular_curves::<2>() {
        for point in curve.points() {
            let order = curve
                .point_order(&point)
                .expect("enumerated point should have a defined order");
            assert_eq!(curve.mul_scalar(&point, order as u64), Ok(curve.identity()));
            assert_eq!(curve.point_has_exact_order(&point, order), Ok(true));
        }
    }
}

#[test]
fn characteristic_three_point_orders_are_defined_for_every_enumerated_point() {
    for curve in exhaustive_nonsingular_curves::<3>() {
        for point in curve.points() {
            let order = curve
                .point_order(&point)
                .expect("enumerated point should have a defined order");
            assert_eq!(curve.mul_scalar(&point, order as u64), Ok(curve.identity()));
            assert_eq!(curve.point_has_exact_order(&point, order), Ok(true));
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_general_weierstrass_group_law_holds_on_enumerated_points_in_characteristic_two(
        (curve, left, right, n, m) in group_law_case::<2>(),
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
        prop_assert_eq!(curve.sub(&left, &right).expect("subtraction should succeed"), curve.add(&left, &curve.neg(&right)).expect("add with inverse should succeed"));
        prop_assert_eq!(curve.double(&left).expect("doubling should succeed"), curve.add(&left, &left).expect("self-addition should succeed"));
        prop_assert_eq!(scalar_sum, split_scalar);
        prop_assert_eq!(curve.mul_scalar_signed(&left, -(n as i64)).expect("signed multiplication should succeed"), curve.mul_scalar(&curve.neg(&left), n).expect("unsigned multiplication should succeed"));
        prop_assert_eq!(curve.neg(&curve.neg(&left)), left);
    }

    #[test]
    fn property_general_weierstrass_group_law_holds_on_enumerated_points_in_characteristic_three(
        (curve, left, right, n, m) in group_law_case::<3>(),
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
        prop_assert_eq!(curve.sub(&left, &right).expect("subtraction should succeed"), curve.add(&left, &curve.neg(&right)).expect("add with inverse should succeed"));
        prop_assert_eq!(curve.double(&left).expect("doubling should succeed"), curve.add(&left, &left).expect("self-addition should succeed"));
        prop_assert_eq!(scalar_sum, split_scalar);
        prop_assert_eq!(curve.mul_scalar_signed(&left, -(n as i64)).expect("signed multiplication should succeed"), curve.mul_scalar(&curve.neg(&left), n).expect("unsigned multiplication should succeed"));
        prop_assert_eq!(curve.neg(&curve.neg(&left)), left);
    }

    #[test]
    fn property_general_weierstrass_group_law_holds_on_enumerated_points_in_characteristic_greater_than_three(
        (curve, left, right, n, m) in group_law_case::<5>(),
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
        prop_assert_eq!(curve.sub(&left, &right).expect("subtraction should succeed"), curve.add(&left, &curve.neg(&right)).expect("add with inverse should succeed"));
        prop_assert_eq!(curve.double(&left).expect("doubling should succeed"), curve.add(&left, &left).expect("self-addition should succeed"));
        prop_assert_eq!(scalar_sum, split_scalar);
        prop_assert_eq!(curve.mul_scalar_signed(&left, -(n as i64)).expect("signed multiplication should succeed"), curve.mul_scalar(&curve.neg(&left), n).expect("unsigned multiplication should succeed"));
        prop_assert_eq!(curve.neg(&curve.neg(&left)), left);
    }

    #[test]
    fn property_general_weierstrass_native_group_law_matches_the_short_companion_when_characteristic_is_supported(
        (curve, left, right, base_scalar, extra_scalar) in group_law_case::<5>(),
        signed_scalar in -7i64..8,
    ) {
        let conversion = curve
            .conversion_to_short_weierstrass()
            .expect("characteristic five should support the short reduction");
        let short_left = conversion
            .map_source_point(&left)
            .expect("sampled point should transport to short");
        let short_right = conversion
            .map_source_point(&right)
            .expect("sampled point should transport to short");
        let scalar = base_scalar + extra_scalar;

        prop_assert_eq!(
            curve.neg(&left),
            conversion
                .map_target_point(&conversion.target().neg(&short_left))
                .expect("negated point should transport back"),
        );
        prop_assert_eq!(
            curve.add(&left, &right).expect("sampled points should add"),
            conversion
                .map_target_point(&conversion.target().add(&short_left, &short_right).expect("short addition should succeed"))
                .expect("sum should transport back"),
        );
        prop_assert_eq!(
            curve.sub(&left, &right).expect("sampled points should subtract"),
            conversion
                .map_target_point(&conversion.target().sub(&short_left, &short_right).expect("short subtraction should succeed"))
                .expect("difference should transport back"),
        );
        prop_assert_eq!(
            curve.double(&left).expect("sampled point should double"),
            conversion
                .map_target_point(&conversion.target().double(&short_left).expect("short doubling should succeed"))
                .expect("double should transport back"),
        );
        prop_assert_eq!(
            curve.mul_scalar(&left, scalar).expect("sampled point should multiply"),
            conversion
                .map_target_point(&conversion.target().mul_scalar(&short_left, scalar).expect("short scalar multiplication should succeed"))
                .expect("multiple should transport back"),
        );
        prop_assert_eq!(
            curve.mul_scalar_signed(&left, signed_scalar).expect("sampled point should multiply by signed scalar"),
            conversion
                .map_target_point(&conversion.target().mul_scalar_signed(&short_left, signed_scalar).expect("short signed scalar multiplication should succeed"))
                .expect("signed multiple should transport back"),
        );
    }
}
