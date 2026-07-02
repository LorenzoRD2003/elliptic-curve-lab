use proptest::prelude::*;

use super::shared::{curve, curve_and_point};
use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::isogenies::{
    error::{IsogenyConstructionError, IsogenyError},
    scalar_multiplication::ScalarMultiplicationIsogeny,
    traits::{Isogeny, VerifiableIsogeny},
};

#[test]
fn constructor_rejects_zero_scalar() {
    assert!(matches!(
        ScalarMultiplicationIsogeny::new(curve(), 0),
        Err(IsogenyError::Construction(
            IsogenyConstructionError::ZeroScalarIsNotIsogeny
        ))
    ));
}

#[test]
fn degree_of_multiplication_by_two_is_four() {
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

    assert_eq!(isogeny.degree(), 4);
    assert_eq!(isogeny.scalar(), &num_bigint::BigUint::from(2u8));
}

#[test]
fn evaluation_matches_group_scalar_multiplication() {
    let curve = curve();
    let point = curve
        .point(
            crate::fields::Fp41::from_i64(3),
            crate::fields::Fp41::from_i64(6),
        )
        .expect("sample point should lie on the curve");
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 3).expect("scalar isogeny should build");

    assert_eq!(
        isogeny.evaluate(&point),
        curve.mul_scalar(&point, 3).map_err(Into::into)
    );
}

#[test]
fn scalar_one_is_identity_map() {
    let curve = curve();
    let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 1)
        .expect("scalar-one isogeny should build");

    for point in curve.points() {
        assert_eq!(
            isogeny
                .evaluate(&point)
                .expect("scalar-one isogeny should evaluate"),
            point
        );
    }
}

#[test]
fn kernel_points_match_the_rational_two_torsion_plus_identity() {
    let curve = curve();
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 2).expect("scalar isogeny should build");

    let mut expected = vec![curve.identity()];
    expected.extend(curve.points_of_order(2));

    assert_eq!(isogeny.kernel_points(), expected.as_slice());
}

#[test]
fn exhaustive_verifier_passes_for_multiplication_by_two() {
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

    assert_eq!(isogeny.verify_maps_domain_to_codomain(), Ok(()));
    assert_eq!(isogeny.verify_maps_kernel_to_identity(), Ok(()));
    assert_eq!(isogeny.verify_homomorphism(), Ok(()));
    assert_eq!(isogeny.verify_kernel_exactness(), Ok(()));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_scalar_isogeny_evaluation_matches_curve_scalar_multiplication(
        (curve, point) in curve_and_point(),
        scalar in 1u64..6,
    ) {
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), scalar)
            .expect("scalar isogeny should build");

        prop_assert_eq!(
            isogeny.evaluate(&point),
            curve.mul_scalar(&point, scalar).map_err(Into::into)
        );
    }

    #[test]
    fn property_scalar_isogeny_kernel_matches_points_killed_by_the_scalar(
        scalar in 1u64..6,
    ) {
        let curve = curve();
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), scalar)
            .expect("scalar isogeny should build");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter(|point| curve.mul_scalar(point, scalar).ok() == Some(curve.identity()))
            .collect();

        prop_assert_eq!(isogeny.kernel_points(), expected.as_slice());
    }
}
