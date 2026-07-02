use crate::fields::traits::*;
use proptest::prelude::*;

use super::shared::F5;
use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    traits::{CurveModel, CurveModelConversion, EnumerableCurveModel},
};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_montgomery_curve,
};

fn reducible_transport_case()
-> impl Strategy<Value = (MontgomeryCurve<F5>, AffinePoint<F5>, AffinePoint<F5>)> {
    arb_nonsingular_montgomery_curve::<crate::fields::Fp5>(CurveStrategyConfig::default())
        .prop_flat_map(|curve| {
            let points = curve.points();
            let len = points.len();

            (Just(curve.clone()), Just(points), 0usize..len, 0usize..len).prop_map(
                |(curve, points, left_index, right_index)| {
                    (
                        curve,
                        points[left_index].clone(),
                        points[right_index].clone(),
                    )
                },
            )
        })
}

#[test]
fn short_reduction_preserves_classical_invariants() {
    let curve = MontgomeryCurve::<F5>::new(F5::one(), F5::one()).expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert!(F5::eq(&curve.c4(), &conversion.target().c4()));
    assert!(F5::eq(&curve.c6(), &conversion.target().c6()));
    assert!(F5::eq(
        &curve.discriminant(),
        &conversion.target().discriminant()
    ));
    assert!(F5::eq(
        &curve.j_invariant(),
        &conversion.target().j_invariant()
    ));
}

#[test]
fn transport_preserves_membership_and_roundtrips_for_all_enumerated_points() {
    let curve = MontgomeryCurve::<F5>::new(F5::one(), F5::one()).expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    for montgomery_point in curve.points() {
        let short_point = conversion
            .map_source_point(&montgomery_point)
            .expect("enumerated Montgomery point should transport to short");

        assert!(curve.contains(&montgomery_point));
        assert!(conversion.target().contains(&short_point));
        assert_eq!(
            conversion
                .map_target_point(&short_point)
                .expect("transported short point should return to the source model"),
            montgomery_point,
        );
    }

    for short_point in conversion.target().points() {
        let montgomery_point = conversion
            .map_target_point(&short_point)
            .expect("enumerated short point should transport to Montgomery");

        assert!(conversion.target().contains(&short_point));
        assert!(curve.contains(&montgomery_point));
        assert_eq!(
            conversion
                .map_source_point(&montgomery_point)
                .expect("transported Montgomery point should return to the target model"),
            short_point,
        );
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_transport_compatibility_matches_the_short_companion_on_random_reducible_curves(
        (curve, left, right) in reducible_transport_case(),
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

        prop_assert!(curve.contains(&left));
        prop_assert!(curve.contains(&right));
        prop_assert!(conversion.target().contains(&short_left));
        prop_assert!(conversion.target().contains(&short_right));
        prop_assert_eq!(curve.j_invariant(), conversion.target().j_invariant());
        prop_assert_eq!(
            conversion
                .map_target_point(&short_left)
                .expect("transported short point should return to Montgomery"),
            left,
        );
        prop_assert_eq!(
            conversion
                .map_target_point(&short_right)
                .expect("transported short point should return to Montgomery"),
            right,
        );
    }
}
