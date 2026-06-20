use crate::elliptic_curves::{
    CurveError, GeneralWeierstrassCurve, ShortWeierstrassCurve,
    traits::{CurveModelConversion, CurveModelConversionError},
};
use crate::fields::{Fp, traits::Field};

type F2 = Fp<2>;
type F5 = Fp<5>;
type F7 = Fp<7>;

#[test]
fn constructor_rejects_singular_coefficients() {
    assert!(matches!(
        GeneralWeierstrassCurve::<F5>::new(
            F5::zero(),
            F5::zero(),
            F5::zero(),
            F5::zero(),
            F5::zero(),
        ),
        Err(CurveError::SingularCurve),
    ));
}

#[test]
fn constructor_allows_characteristic_two_when_discriminant_is_nonzero() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular general Weierstrass curve in characteristic two");

    assert!(F2::eq(curve.a1(), &F2::one()));
    assert!(F2::eq(&curve.discriminant(), &F2::one()));
}

#[test]
fn short_embedding_example_matches_the_expected_classical_invariants() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    assert!(F7::eq(curve.a1(), &F7::zero()));
    assert!(F7::eq(curve.a2(), &F7::zero()));
    assert!(F7::eq(curve.a3(), &F7::zero()));
    assert!(F7::eq(curve.a4(), &F7::from_i64(2)));
    assert!(F7::eq(curve.a6(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.b2(), &F7::zero()));
    assert!(F7::eq(&curve.b4(), &F7::from_i64(4)));
    assert!(F7::eq(&curve.b6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.b8(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c4(), &F7::from_i64(2)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(5)));
}

#[test]
fn weierstrass_invariants_satisfy_the_classical_relation() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
    let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

    assert!(F7::eq(&left, &right));
}

#[test]
fn display_and_debug_surface_the_general_equation() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let display = curve.to_string();

    assert!(display.starts_with("y^2 + ("));
    assert!(display.contains(")xy + ("));
    assert!(display.contains(")y = x^3 + ("));
    assert!(display.contains(")x^2 + ("));
    assert!(display.contains(")x + ("));
    assert!(format!("{curve:?}").contains("GeneralWeierstrassCurve"));
}

#[test]
fn clone_and_equality_preserve_all_coefficients() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let clone = curve.clone();

    assert_eq!(clone, curve);
    assert!(F7::eq(clone.a1(), curve.a1()));
    assert!(F7::eq(clone.a2(), curve.a2()));
    assert!(F7::eq(clone.a3(), curve.a3()));
    assert!(F7::eq(clone.a4(), curve.a4()));
    assert!(F7::eq(clone.a6(), curve.a6()));
}

#[test]
fn equation_string_mentions_every_general_weierstrass_term() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let equation = curve.to_equation_string();

    assert!(equation.contains("xy"));
    assert!(equation.contains("x^2"));
    assert!(equation.contains("x^3"));
    assert!(equation.contains("y = x^3"));
}

#[test]
fn reduction_rejects_characteristic_three() {
    type F3 = Fp<3>;

    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("the general model itself may exist in characteristic three");

    assert!(matches!(
        curve.conversion_to_short_weierstrass(),
        Err(CurveModelConversionError::UnsupportedCharacteristic { characteristic: 3 }),
    ));
}

#[test]
fn reduction_rejects_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("the general model itself may exist in characteristic two");

    assert!(matches!(
        curve.conversion_to_short_weierstrass(),
        Err(CurveModelConversionError::UnsupportedCharacteristic { characteristic: 2 }),
    ));
}

#[test]
fn reduction_of_an_embedded_short_curve_has_zero_transport_parameters() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic seven should support the reduction");

    assert_eq!(conversion.source(), &curve);
    assert!(F7::eq(conversion.target().a(), &F7::from_i64(2)));
    assert!(F7::eq(conversion.target().b(), &F7::from_i64(3)));
}

#[test]
fn reduction_produces_expected_short_companion_for_a_genuinely_general_example() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    assert_eq!(conversion.source(), &curve);
    assert!(F5::eq(conversion.target().a(), &F5::from_i64(4)));
    assert!(F5::eq(conversion.target().b(), &F5::from_i64(4)));
}

#[test]
fn try_as_short_weierstrass_matches_the_reduction_companion() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let from_helper = curve
        .try_as_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let from_reduction = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction")
        .target()
        .clone();

    assert_eq!(from_helper, from_reduction);
}

#[test]
fn try_from_general_weierstrass_reference_matches_the_reduction_companion() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let from_try_from = ShortWeierstrassCurve::try_from(&curve)
        .expect("characteristic five should support the reduction");
    let from_reduction = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction")
        .target()
        .clone();

    assert_eq!(from_try_from, from_reduction);
}

#[test]
fn transporting_infinity_between_general_and_short_models_is_stable() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let short_infinity = conversion
        .map_source_point(&crate::elliptic_curves::AffinePoint::<F5>::Infinity)
        .expect("infinity should transport to the short companion");
    let general_infinity = conversion
        .map_target_point(&crate::elliptic_curves::AffinePoint::<F5>::Infinity)
        .expect("infinity should transport back to the general model");

    assert_eq!(
        short_infinity,
        crate::elliptic_curves::AffinePoint::<F5>::Infinity
    );
    assert_eq!(
        general_infinity,
        crate::elliptic_curves::AffinePoint::<F5>::Infinity
    );
}

#[test]
fn transporting_a_general_point_to_short_matches_the_expected_coordinates() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let general_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero());

    let image = conversion
        .map_source_point(&general_point)
        .expect("point should transport to the short companion");

    assert_eq!(
        image,
        crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::from_i64(3))
    );
}

#[test]
fn transporting_short_and_general_points_roundtrips() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let general_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero());
    let short_point = conversion
        .map_source_point(&general_point)
        .expect("point should transport to the short companion");

    let general_roundtrip = conversion
        .map_target_point(&short_point)
        .expect("short point should transport back to the general model");
    let short_roundtrip = conversion
        .map_source_point(&general_roundtrip)
        .expect("transporting back again should still succeed");

    assert_eq!(general_roundtrip, general_point);
    assert_eq!(short_roundtrip, short_point);
}

#[test]
fn short_reduction_to_general_embeds_the_short_coefficients() {
    let short_curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let conversion = short_curve.conversion_to_general_weierstrass();

    assert_eq!(conversion.source(), &short_curve);
    assert!(F7::eq(conversion.target().a1(), &F7::zero()));
    assert!(F7::eq(conversion.target().a2(), &F7::zero()));
    assert!(F7::eq(conversion.target().a3(), &F7::zero()));
    assert!(F7::eq(conversion.target().a4(), &F7::from_i64(2)));
    assert!(F7::eq(conversion.target().a6(), &F7::from_i64(3)));
}

#[test]
fn as_general_weierstrass_matches_the_reduction_general_curve() {
    let short_curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let from_helper = short_curve.as_general_weierstrass();
    let from_reduction = short_curve
        .conversion_to_general_weierstrass()
        .target()
        .clone();

    assert_eq!(from_helper, from_reduction);
}

#[test]
fn from_short_weierstrass_reference_matches_the_reduction_general_curve() {
    let short_curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    let from_trait = GeneralWeierstrassCurve::from(&short_curve);
    let from_reduction = short_curve
        .conversion_to_general_weierstrass()
        .target()
        .clone();

    assert_eq!(from_trait, from_reduction);
}

#[test]
fn transporting_a_short_point_to_general_is_the_identity_embedding() {
    let short_curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let conversion = short_curve.conversion_to_general_weierstrass();
    let short_point =
        crate::elliptic_curves::AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(1));

    let general_point = conversion
        .map_source_point(&short_point)
        .expect("embedded short point should lie on the general model");
    let short_roundtrip = conversion
        .map_target_point(&general_point)
        .expect("embedded general point should return to the short model");

    assert_eq!(general_point, short_point);
    assert_eq!(short_roundtrip, short_point);
}

#[test]
fn conversion_reports_invalid_source_and_target_points_honestly() {
    let general_curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = general_curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the conversion");
    let bad_general_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::one());
    let bad_short_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::one(), F5::one());

    assert_eq!(
        conversion.map_source_point(&bad_general_point),
        Err(CurveModelConversionError::PointNotOnSource)
    );
    assert_eq!(
        conversion.map_target_point(&bad_short_point),
        Err(CurveModelConversionError::PointNotOnTarget)
    );
}
