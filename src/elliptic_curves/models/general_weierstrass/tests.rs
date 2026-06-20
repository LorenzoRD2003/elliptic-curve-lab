use num_complex::Complex64;

use crate::elliptic_curves::{
    CurveError, GeneralWeierstrassCurve, ShortWeierstrassCurve,
    traits::{
        AffineCurveModel, CurveModel, CurveModelConversion, CurveModelConversionError,
        FiniteGroupCurveModel, GroupCurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints,
    },
};
use crate::fields::{
    ComplexApprox, Fp, Q,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
    traits::Field,
};

type F2 = Fp<2>;
type F5 = Fp<5>;
type F7 = Fp<7>;

#[derive(Clone, Copy)]
struct F4GeneralWeierstrassSpec;

impl ExtensionFieldSpec for F4GeneralWeierstrassSpec {
    type Base = F2;

    const NAME: &'static str = "F4 for general Weierstrass lifting tests";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![F2::one(), F2::one(), F2::one()])
            .expect("x^2 + x + 1 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

type F4 = ExtensionField<F4GeneralWeierstrassSpec>;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

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
fn curve_model_identity_and_membership_helpers_work_for_general_weierstrass() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let finite_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero());
    let infinity = crate::elliptic_curves::AffinePoint::<F5>::infinity();
    let off_curve_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::one());

    assert_eq!(curve.identity(), infinity);
    assert!(curve.contains(&finite_point));
    assert!(curve.is_on_curve_nonzero(&finite_point));
    assert!(curve.contains(&infinity));
    assert!(curve.is_identity(&infinity));
    assert!(!curve.is_on_curve_nonzero(&infinity));
    assert!(!curve.contains(&off_curve_point));
    assert!(!curve.is_on_curve_nonzero(&off_curve_point));
}

#[test]
fn affine_curve_model_point_accepts_valid_general_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let point = curve
        .point(F5::zero(), F5::zero())
        .expect("the affine point should lie on the curve");

    assert_eq!(
        point,
        crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero())
    );
    assert!(curve.contains(&point));
}

#[test]
fn affine_curve_model_point_rejects_invalid_general_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    assert_eq!(
        curve.point(F5::zero(), F5::one()),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_general_weierstrass_invariant() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    assert!(F7::eq(
        &HasJInvariant::j_invariant(&curve),
        &GeneralWeierstrassCurve::j_invariant(&curve)
    ));
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
fn affine_curve_model_point_accepts_valid_general_points_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    let point = curve
        .point(F2::one(), F2::zero())
        .expect("the affine point should lie on the characteristic-two curve");

    assert!(curve.contains(&point));
}

#[test]
fn general_weierstrass_negation_uses_the_model_specific_involution() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::zero(), F5::from_i64(4))
        .expect("sample point should lie on the curve");

    assert_eq!(
        curve.neg(&point),
        curve
            .point(F5::zero(), F5::zero())
            .expect("the inverse should lie on the curve")
    );
    assert_ne!(curve.neg(&point), point.neg());
}

#[test]
fn general_weierstrass_addition_matches_the_short_companion_transport() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let left = curve
        .point(F5::zero(), F5::zero())
        .expect("left point should lie on the curve");
    let right = curve
        .point(F5::from_i64(2), F5::one())
        .expect("right point should lie on the curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let expected = conversion
        .map_target_point(
            &conversion
                .target()
                .add(
                    &conversion
                        .map_source_point(&left)
                        .expect("left point should transport to short"),
                    &conversion
                        .map_source_point(&right)
                        .expect("right point should transport to short"),
                )
                .expect("short companion addition should succeed"),
        )
        .expect("short sum should transport back");

    assert_eq!(curve.add(&left, &right), Ok(expected));
}

#[test]
fn general_weierstrass_doubling_matches_the_short_companion_transport() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::zero(), F5::zero())
        .expect("sample point should lie on the curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let expected = conversion
        .map_target_point(
            &conversion
                .target()
                .double(
                    &conversion
                        .map_source_point(&point)
                        .expect("point should transport to short"),
                )
                .expect("short companion doubling should succeed"),
        )
        .expect("short double should transport back");

    assert_eq!(curve.double(&point), Ok(expected));
}

#[test]
fn general_weierstrass_small_torsion_helpers_work_with_the_native_affine_group_law() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::from_i64(2), F5::one())
        .expect("sample point should lie on the curve");

    assert_eq!(curve.double(&point), Ok(curve.identity()));
    assert_eq!(curve.point_has_exact_order(&point, 2), Ok(true));
    assert!(curve.is_torsion_point(&point, 2));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(curve.identity()));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn general_weierstrass_group_law_handles_characteristic_two_natively() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let point = curve
        .point(F2::one(), F2::zero())
        .expect("sample point should lie on the curve");

    assert_eq!(curve.add(&point, &point), Ok(curve.identity()));
    assert_eq!(curve.neg(&point), point);
    assert_eq!(curve.double(&point), Ok(curve.identity()));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(curve.identity()));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn general_weierstrass_group_law_handles_characteristic_three_natively() {
    type F3 = Fp<3>;

    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let point = curve
        .point(F3::zero(), F3::zero())
        .expect("sample point should lie on the curve");
    let inverse = curve
        .point(F3::zero(), F3::from_i64(2))
        .expect("the inverse point should lie on the curve");
    let doubled = curve
        .point(F3::from_i64(2), F3::one())
        .expect("the doubled point should lie on the curve");

    assert_eq!(curve.neg(&point), inverse);
    assert_eq!(curve.add(&point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.double(&point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(doubled));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn lift_x_over_an_odd_characteristic_curve_returns_two_points_when_the_fiber_is_split() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let lifted = curve
        .lift_x(F5::zero())
        .expect("finite-field odd-characteristic lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points, got {other:?}"),
    }
}

#[test]
fn lift_x_over_an_odd_characteristic_curve_returns_no_point_when_the_fiber_is_empty() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    assert_eq!(
        curve
            .lift_x(F5::from_i64(3))
            .expect("finite-field odd-characteristic lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_one_point_when_b_is_zero() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    assert_eq!(
        curve
            .lift_x(F2::one())
            .expect("characteristic-two lifting should succeed"),
        LiftedPoints::OnePoint(crate::elliptic_curves::AffinePoint::<F2>::new(
            F2::one(),
            F2::zero()
        ))
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_no_point_when_the_artin_schreier_equation_is_unsolvable() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    assert_eq!(
        curve
            .lift_x(F2::zero())
            .expect("characteristic-two lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_two_points_when_the_artin_schreier_equation_is_solvable() {
    let curve = GeneralWeierstrassCurve::<F4>::new(
        F4::zero(),
        F4::zero(),
        F4::one(),
        F4::zero(),
        F4::zero(),
    )
    .expect("non-singular curve in characteristic two");

    let lifted = curve
        .lift_x(F4::one())
        .expect("characteristic-two lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points, got {other:?}"),
    }
}

#[test]
fn lift_x_over_q_returns_two_rational_points_when_the_fiber_is_split() {
    let curve =
        GeneralWeierstrassCurve::<Q>::new(Q::one(), Q::one(), Q::one(), Q::one(), Q::zero())
            .expect("non-singular curve over Q");

    let lifted = curve
        .lift_x(Q::zero())
        .expect("rational lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points over Q, got {other:?}"),
    }
}

#[test]
fn lift_x_over_q_returns_no_point_when_the_quadratic_has_no_rational_root() {
    let curve =
        GeneralWeierstrassCurve::<Q>::new(Q::one(), Q::one(), Q::one(), Q::one(), Q::zero())
            .expect("non-singular curve over Q");

    assert_eq!(
        curve
            .lift_x(Q::from_i64(3))
            .expect("rational lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_over_complex_approx_returns_two_points() {
    let curve = GeneralWeierstrassCurve::<ComplexApprox>::new(
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(0.0, 0.0),
    )
    .expect("non-singular curve over ComplexApprox");

    let lifted = curve
        .lift_x(c(3.0, 0.0))
        .expect("complex lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points over ComplexApprox, got {other:?}"),
    }
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
