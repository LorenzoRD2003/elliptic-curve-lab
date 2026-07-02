use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    traits::{
        CurveModel, CurveModelConversion, EnumerableCurveModel, FiniteGroupCurveModel,
        GroupCurveModel,
    },
};
use crate::fields::traits::EnumerableFiniteField;
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_montgomery_curve,
};

fn exhaustive_nonsingular_curves<F>() -> Vec<MontgomeryCurve<F>>
where
    F: EnumerableFiniteField,
{
    let elements = F::elements();
    let mut curves = Vec::new();
    for a in &elements {
        for b in &elements {
            if let Ok(curve) = MontgomeryCurve::<F>::new(a.clone(), b.clone()) {
                curves.push(curve);
            }
        }
    }
    curves
}

fn group_law_case<F>()
-> impl Strategy<Value = (MontgomeryCurve<F>, AffinePoint<F>, AffinePoint<F>, u64, u64)>
where
    F: EnumerableFiniteField + crate::fields::traits::SqrtField + 'static,
    F::Elem: 'static,
{
    arb_nonsingular_montgomery_curve::<F>(CurveStrategyConfig::default()).prop_flat_map(|curve| {
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

#[test]
fn exhaustive_group_axiom_check_passes_for_every_nonsingular_curve_over_f3() {
    for curve in exhaustive_nonsingular_curves::<crate::fields::Fp3>() {
        assert_eq!(curve.check_group_axioms(), Ok(()));
    }
}

#[test]
fn exhaustive_group_axiom_check_passes_for_every_nonsingular_curve_over_f5() {
    for curve in exhaustive_nonsingular_curves::<crate::fields::Fp5>() {
        assert_eq!(curve.check_group_axioms(), Ok(()));
    }
}

#[test]
fn characteristic_three_point_orders_are_defined_for_every_enumerated_point() {
    for curve in exhaustive_nonsingular_curves::<crate::fields::Fp3>() {
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
fn characteristic_five_point_orders_are_defined_for_every_enumerated_point() {
    for curve in exhaustive_nonsingular_curves::<crate::fields::Fp5>() {
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
    fn property_montgomery_group_law_holds_on_enumerated_points_in_characteristic_three(
        (curve, left, right, n, m) in group_law_case::<crate::fields::Fp3>(),
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
    fn property_montgomery_group_law_holds_on_enumerated_points_in_characteristic_five(
        (curve, left, right, n, m) in group_law_case::<crate::fields::Fp5>(),
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
    fn property_montgomery_native_group_law_matches_the_short_companion_when_characteristic_is_supported(
        (curve, left, right, base_scalar, extra_scalar) in group_law_case::<crate::fields::Fp5>(),
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
