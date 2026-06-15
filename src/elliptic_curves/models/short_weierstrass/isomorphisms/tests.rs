use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::{
        CurveIsomorphismError, ShortWeierstrassIsomorphism, ShortWeierstrassQuadraticTwist,
        TwistKind,
    },
    traits::{AffineCurveModel, CurveIsomorphism, CurveModel},
};
use crate::fields::{
    FieldError, Fp,
    extension_field::define_fp_quadratic_extension,
    traits::{CbrtField, EnumerableFiniteField, Field, SqrtField},
};
use crate::proptest_support::{
    config::CurveStrategyConfig, isogenies::arb_short_weierstrass_isomorphism_case,
};

type F7 = Fp<7>;
type F13 = Fp<13>;
type F19 = Fp<19>;

define_fp_quadratic_extension!(
    spec: F7Sqrt3Spec,
    field: F7Sqrt3,
    base: F7,
    non_residue: 3,
    name: "F7(sqrt(3))",
);

#[allow(dead_code)]
type _F7Sqrt3Marker = F7Sqrt3;

define_fp_quadratic_extension!(
    spec: F19Sqrt2Spec,
    field: F19Sqrt2,
    base: F19,
    non_residue: 2,
    name: "F19(sqrt(2))",
);

define_fp_quadratic_extension!(
    spec: F19Sqrt3Spec,
    field: F19Sqrt3,
    base: F19,
    non_residue: 3,
    name: "F19(sqrt(3))",
);

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

fn f19_curve() -> ShortWeierstrassCurve<F19> {
    ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
}

fn f7_j1728_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::zero()).expect("valid j=1728 curve")
}

fn f13_j0_curve() -> ShortWeierstrassCurve<F13> {
    ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(1)).expect("valid j=0 curve")
}

fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    f7_curve()
        .point(F7::from_i64(x), F7::from_i64(y))
        .expect("valid point on the sample curve")
}

fn f19_point(x: i64, y: i64) -> AffinePoint<F19> {
    f19_curve()
        .point(F19::from_i64(x), F19::from_i64(y))
        .expect("valid point on the sample curve")
}

#[test]
fn constructor_rejects_noninvertible_scale() {
    assert!(matches!(
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::zero()),
        Err(CurveIsomorphismError::NonInvertibleScale)
    ));
}

#[test]
fn domain_getter_returns_original_curve() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

    assert!(F7::eq(isomorphism.domain().a(), &F7::from_i64(2)));
    assert!(F7::eq(isomorphism.domain().b(), &F7::from_i64(3)));
}

#[test]
fn codomain_is_derived_from_u4_and_u6() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let codomain = isomorphism.codomain();

    assert!(F7::eq(codomain.a(), &F7::from_i64(1)));
    assert!(F7::eq(codomain.b(), &F7::from_i64(3)));
}

#[test]
fn scaling_factor_getter_returns_u() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

    assert!(F7::eq(isomorphism.scaling_factor(), &F7::from_i64(3)));
}

#[test]
fn inverse_uses_inverse_scaling_factor() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let inverse = isomorphism.inverse().expect("inverse should exist");

    assert!(F7::eq(inverse.scaling_factor(), &F7::from_i64(5)));
}

#[test]
fn inverse_domain_is_the_original_codomain() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let codomain = isomorphism.codomain();
    let inverse = isomorphism.inverse().expect("inverse should exist");

    assert!(F7::eq(inverse.domain().a(), codomain.a()));
    assert!(F7::eq(inverse.domain().b(), codomain.b()));
}

#[test]
fn inverse_codomain_is_the_original_domain() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let inverse = isomorphism.inverse().expect("inverse should exist");
    let inverse_codomain = inverse.codomain();

    assert!(F7::eq(inverse_codomain.a(), isomorphism.domain().a()));
    assert!(F7::eq(inverse_codomain.b(), isomorphism.domain().b()));
}

#[test]
fn inverse_of_inverse_returns_original_scaling_and_curves() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let double_inverse = isomorphism
        .inverse()
        .expect("inverse should exist")
        .inverse()
        .expect("double inverse should exist");

    assert!(F7::eq(
        double_inverse.scaling_factor(),
        isomorphism.scaling_factor()
    ));
    assert!(F7::eq(
        double_inverse.domain().a(),
        isomorphism.domain().a()
    ));
    assert!(F7::eq(
        double_inverse.domain().b(),
        isomorphism.domain().b()
    ));

    let double_inverse_codomain = double_inverse.codomain();
    let original_codomain = isomorphism.codomain();
    assert!(F7::eq(double_inverse_codomain.a(), original_codomain.a()));
    assert!(F7::eq(double_inverse_codomain.b(), original_codomain.b()));
}

#[test]
fn evaluate_sends_infinity_to_infinity() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

    assert_eq!(
        isomorphism
            .evaluate(&AffinePoint::infinity())
            .expect("the identity should map to itself"),
        AffinePoint::infinity()
    );
}

#[test]
fn evaluate_rejects_point_outside_the_domain() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert!(matches!(
        isomorphism.evaluate(&invalid),
        Err(CurveIsomorphismError::PointNotOnDomain)
    ));
}

#[test]
fn evaluate_transports_finite_points_to_the_codomain() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let point = f7_point(2, 1);
    let image = isomorphism
        .evaluate(&point)
        .expect("a domain point should map into the codomain");

    assert!(isomorphism.codomain().contains(&image));
    assert_eq!(image, AffinePoint::new(F7::from_i64(4), F7::from_i64(6)));
}

#[test]
fn inverse_recovers_the_original_point_after_evaluation() {
    let isomorphism =
        ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
    let point = f7_point(2, 1);
    let image = isomorphism
        .evaluate(&point)
        .expect("a domain point should map into the codomain");

    assert_eq!(
        isomorphism
            .inverse()
            .expect("inverse should exist")
            .evaluate(&image)
            .expect("the inverse should recover the original point"),
        point
    );
}

#[test]
fn quadratic_twist_package_stores_original_twist_and_factor() {
    let original = f19_curve();
    let package = ShortWeierstrassQuadraticTwist::new(original, F19::from_i64(2))
        .expect("non-zero twist factor should produce a valid package");

    assert!(F19::eq(package.original().a(), &F19::from_i64(2)));
    assert!(F19::eq(package.original().b(), &F19::from_i64(3)));
    assert!(F19::eq(package.factor(), &F19::from_i64(2)));
    assert!(package.original().has_same_j_invariant(package.twist()));
}

#[test]
fn quadratic_twist_kind_is_trivial_when_the_factor_is_a_square() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
        .expect("square twist factor should produce a valid package");

    assert_eq!(package.kind(), TwistKind::Trivial);
}

#[test]
fn quadratic_twist_kind_is_quadratic_when_the_factor_is_not_a_square() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
        .expect("non-square twist factor should produce a valid package");

    assert_eq!(package.kind(), TwistKind::Quadratic);
}

#[test]
fn base_field_isomorphism_exists_for_generic_square_and_non_square_factors() {
    let trivial_package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
        .expect("square twist factor should produce a valid package");
    let quadratic_package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
        .expect("non-square twist factor should produce a valid package");

    assert!(trivial_package.base_field_isomorphism().is_some());
    assert!(quadratic_package.base_field_isomorphism().is_none());
}

#[test]
fn j1728_non_square_factor_can_still_be_trivial() {
    let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
        .expect("non-square twist factor should still produce a valid package");

    assert_eq!(package.kind(), TwistKind::Trivial);
}

#[test]
fn j1728_non_square_factor_can_still_produce_a_base_field_isomorphism() {
    let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
        .expect("non-square twist factor should still produce a valid package");
    let isomorphism = package
        .base_field_isomorphism()
        .expect("j = 1728 should admit the extra base-field witness here");

    assert!(F7::eq(isomorphism.scaling_factor(), &F7::from_i64(2)));
    assert!(F7::eq(isomorphism.codomain().a(), package.twist().a()));
    assert!(F7::eq(isomorphism.codomain().b(), package.twist().b()));
}

#[test]
fn quadratic_extension_isomorphism_rejects_j1728_twists_that_are_already_trivial() {
    let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
        .expect("non-square twist factor should still produce a valid package");

    assert!(matches!(
        package.isomorphism_over_quadratic_extension::<F7Sqrt3Spec>(),
        Err(CurveIsomorphismError::Field(
            FieldError::NonIrreduciblePolynomial
        ))
    ));
}

#[test]
fn j0_square_factor_is_still_trivial() {
    let package = ShortWeierstrassQuadraticTwist::new(f13_j0_curve(), F13::from_i64(4))
        .expect("square twist factor should produce a valid package");

    assert_eq!(package.kind(), TwistKind::Trivial);
    assert!(package.base_field_isomorphism().is_some());
}

#[test]
fn j0_non_square_factor_stays_quadratic_in_the_sample_prime_field() {
    let package = ShortWeierstrassQuadraticTwist::new(f13_j0_curve(), F13::from_i64(2))
        .expect("non-square twist factor should produce a valid package");

    assert_eq!(package.kind(), TwistKind::Quadratic);
    assert!(package.base_field_isomorphism().is_none());
}

#[test]
fn j0_sample_field_has_nontrivial_cube_roots_of_unity_but_they_are_still_squares() {
    let nontrivial_cube_roots_of_unity = F13::elements()
        .into_iter()
        .filter(|element| {
            F13::eq(&F13::cube(element), &F13::one()) && !F13::eq(element, &F13::one())
        })
        .collect::<Vec<_>>();

    assert!(!nontrivial_cube_roots_of_unity.is_empty());
    assert!(
        nontrivial_cube_roots_of_unity
            .iter()
            .all(F13::has_square_root)
    );
    assert!(
        nontrivial_cube_roots_of_unity
            .iter()
            .all(|element| !F13::has_cube_root(element))
    );
}

#[test]
fn base_field_isomorphism_transports_points_into_the_stored_twist() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
        .expect("square twist factor should produce a valid package");
    let point = f19_point(3, 6);
    let isomorphism = package
        .base_field_isomorphism()
        .expect("a square twist factor should produce a base-field isomorphism");
    let image = isomorphism
        .evaluate(&point)
        .expect("the base-field isomorphism should transport points");

    assert!(package.twist().contains(&image));
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
fn quadratic_extension_isomorphism_uses_the_expected_extension_and_codomain() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
        .expect("non-square twist factor should produce a valid package");
    let isomorphism = package
        .isomorphism_over_quadratic_extension::<F19Sqrt2Spec>()
        .expect("the matching quadratic extension should produce an isomorphism");
    let lifted_twist = ShortWeierstrassCurve::<F19Sqrt2>::new(
        F19Sqrt2::from_base(*package.twist().a()),
        F19Sqrt2::from_base(*package.twist().b()),
    )
    .expect("the stored twist should lift to the extension field");
    let derived_codomain = isomorphism.codomain();

    assert!(F19Sqrt2::eq(
        &F19Sqrt2::square(isomorphism.scaling_factor()),
        &F19Sqrt2::from_base(F19::from_i64(2))
    ));
    assert!(F19Sqrt2::eq(derived_codomain.a(), lifted_twist.a()));
    assert!(F19Sqrt2::eq(derived_codomain.b(), lifted_twist.b()));
}

#[test]
fn quadratic_extension_isomorphism_transports_points_and_inverse_recovers_them() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
        .expect("non-square twist factor should produce a valid package");
    let isomorphism = package
        .isomorphism_over_quadratic_extension::<F19Sqrt2Spec>()
        .expect("the matching quadratic extension should produce an isomorphism");
    let lifted_domain = ShortWeierstrassCurve::<F19Sqrt2>::new(
        F19Sqrt2::from_base(*package.original().a()),
        F19Sqrt2::from_base(*package.original().b()),
    )
    .expect("the original curve should lift to the extension field");
    let point = lifted_domain
        .point(
            F19Sqrt2::from_base(F19::from_i64(3)),
            F19Sqrt2::from_base(F19::from_i64(6)),
        )
        .expect("the sample point should stay on the lifted curve");
    let image = isomorphism
        .evaluate(&point)
        .expect("the quadratic-extension isomorphism should transport points");

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
fn quadratic_extension_isomorphism_rejects_square_twist_factors() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
        .expect("square twist factor should produce a valid package");

    assert!(matches!(
        package.isomorphism_over_quadratic_extension::<F19Sqrt2Spec>(),
        Err(CurveIsomorphismError::Field(
            FieldError::NonIrreduciblePolynomial
        ))
    ));
}

#[test]
fn quadratic_extension_isomorphism_rejects_incompatible_extension_specs() {
    let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
        .expect("non-square twist factor should produce a valid package");
    let _field = F19Sqrt3::new().expect("the alternative quadratic extension should validate");

    assert!(matches!(
        package.isomorphism_over_quadratic_extension::<F19Sqrt3Spec>(),
        Err(CurveIsomorphismError::Field(
            FieldError::IncompatibleFieldParameters
        ))
    ));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn generated_short_weierstrass_isomorphisms_preserve_domain_points(
        case in arb_short_weierstrass_isomorphism_case::<19>(CurveStrategyConfig::default()),
    ) {
        let image = case
            .isomorphism
            .evaluate(&case.sample_point)
            .expect("generated domain point should evaluate");

        prop_assert!(case.curve.contains(&case.sample_point));
        prop_assert!(case.isomorphism.codomain().contains(&image));
        prop_assert_eq!(
            case.isomorphism
                .inverse()
                .expect("generated isomorphism should be invertible")
                .evaluate(&image)
                .expect("inverse should recover the sampled point"),
            case.sample_point
        );
    }
}
