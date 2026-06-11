use proptest::prelude::*;

use crate::elliptic_curves::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
    ShortWeierstrassCurve, ShortWeierstrassFunction, ShortWeierstrassFunctionField,
};
use crate::fields::{Field, Fp};
use crate::isogenies::scalar_multiplication::ScalarMultiplicationIsogeny;
use crate::isogenies::{
    AbsoluteFrobeniusIsogeny, DualIsogenyError, FrobeniusLikeIsogeny, Isogeny,
    IsogenyConstructionError, IsogenyError, KernelDescription, VerifiableIsogeny,
    VerschiebungCertificate, VerschiebungIsogeny,
};
use crate::polynomials::evaluation::evaluate_dense;

type F41 = Fp<41>;
type Curve = ShortWeierstrassCurve<F41>;

crate::fields::define_fp_quadratic_extension!(
    spec: F5Sqrt2ScalarMultiplicationSpec,
    field: F5Sqrt2ScalarMultiplication,
    base: Fp<5>,
    non_residue: 2,
    name: "F5(sqrt(2)) for scalar-multiplication Frobenius tests",
);

fn curve() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn nontrivial_extension_curve() -> ShortWeierstrassCurve<F5Sqrt2ScalarMultiplication> {
    let alpha = F5Sqrt2ScalarMultiplication::element(vec![Fp::<5>::zero(), Fp::<5>::one()]);
    ShortWeierstrassCurve::<F5Sqrt2ScalarMultiplication>::new(
        alpha,
        F5Sqrt2ScalarMultiplication::one(),
    )
    .expect("valid curve over F5^2")
}

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
    assert_eq!(isogeny.scalar(), 2);
}

#[test]
fn evaluation_matches_group_scalar_multiplication() {
    let curve = curve();
    let point = curve
        .point(F41::from_i64(3), F41::from_i64(6))
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

#[test]
fn function_field_map_from_verschiebung_recovers_the_certified_p_pullback() {
    let curve = curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let candidate_v = crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
        frobenius.codomain().clone(),
        frobenius.domain().clone(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .x(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .y(),
    )
    .expect("identity candidate on the twist should validate");
    let verschiebung =
        VerschiebungIsogeny::new(frobenius.clone(), candidate_v).expect("candidate should build");
    let expected_left = frobenius
        .as_function_field_map()
        .compose(verschiebung.as_function_field_map())
        .expect("left composition should build");
    let expected_right = verschiebung
        .as_function_field_map()
        .compose(&frobenius.as_function_field_map())
        .expect("right composition should build");
    let certificate =
        VerschiebungCertificate::new(verschiebung, expected_left.clone(), expected_right)
            .expect("certificate should build");

    let scalar =
        ScalarMultiplicationIsogeny::new(curve, 41).expect("scalar multiplication should build");

    assert_eq!(
        scalar
            .as_function_field_map_from_verschiebung(&certificate)
            .expect("certified map should build"),
        expected_left
    );
}

#[test]
fn direct_p_pullback_can_build_a_verschiebung_isogeny() {
    let scalar =
        ScalarMultiplicationIsogeny::new(curve(), 41).expect("scalar multiplication should build");
    let verschiebung = scalar
        .verschiebung_isogeny_from_direct_p_pullback()
        .expect("Verschiebung should be extracted from [p]^*");

    assert_eq!(verschiebung.codomain_curve(), scalar.domain());
    assert_eq!(
        verschiebung.domain_curve(),
        verschiebung.frobenius().codomain()
    );
    assert_eq!(verschiebung.degree(), 41);
}

#[test]
fn direct_p_pullback_can_build_a_verified_verschiebung_over_nontrivial_extension_curve() {
    let curve = nontrivial_extension_curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 5)
        .expect("scalar multiplication should build");
    let verschiebung = scalar
        .verschiebung_isogeny_from_direct_p_pullback()
        .expect("Verschiebung should be extracted from [p]^*");
    let certificate = scalar
        .verschiebung_certificate_from_direct_p_pullback()
        .expect("certificate should build");

    assert_eq!(verschiebung.codomain_curve(), &curve);
    assert_eq!(
        verschiebung.domain_curve(),
        verschiebung.frobenius().codomain()
    );
    assert_ne!(verschiebung.domain_curve(), verschiebung.codomain_curve());
    assert_eq!(certificate.verify_duality_relations(), Ok(()));
}

#[test]
fn function_field_map_of_scalar_one_is_the_identity_pullback() {
    let curve = curve();
    let field = ShortWeierstrassFunctionField::<F41>::new(curve.clone());
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 1)
        .expect("scalar multiplication should build");

    assert_eq!(
        scalar
            .as_function_field_map()
            .expect("function-field pullback should build"),
        crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
            curve.clone(),
            curve,
            field.x(),
            field.y(),
        )
        .expect("identity pullback should validate")
    );
}

#[test]
fn direct_function_field_map_matches_point_evaluation_away_from_the_kernel() {
    let curve = curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 2)
        .expect("scalar multiplication should build");
    let map = scalar
        .as_function_field_map()
        .expect("function-field pullback should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !scalar.kernel_points().contains(point))
        .expect("sample curve should have a point outside the kernel");
    let image = scalar
        .evaluate(&point)
        .expect("sample point should evaluate under [2]");

    let x_value = evaluate_short_weierstrass_function_at_point(map.x_pullback(), &point)
        .expect("non-kernel point should avoid poles in x pullback");
    let y_value = evaluate_short_weierstrass_function_at_point(map.y_pullback(), &point)
        .expect("non-kernel point should avoid poles in y pullback");

    assert_eq!(
        Some(&x_value),
        crate::elliptic_curves::AffinePoint::x_coordinate(&image)
    );
    assert_eq!(
        Some(&y_value),
        crate::elliptic_curves::AffinePoint::y_coordinate(&image)
    );
}

#[test]
fn direct_p_pullback_matches_point_evaluation_away_from_the_kernel() {
    let curve = curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 41)
        .expect("scalar multiplication should build");
    let map = scalar
        .as_function_field_map()
        .expect("direct [p]^* should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !scalar.kernel_points().contains(point))
        .expect("sample curve should have a point outside the kernel");
    let image = scalar
        .evaluate(&point)
        .expect("sample point should evaluate under [p]");

    let x_value = evaluate_short_weierstrass_function_at_point(map.x_pullback(), &point)
        .expect("non-kernel point should avoid poles in x pullback");
    let y_value = evaluate_short_weierstrass_function_at_point(map.y_pullback(), &point)
        .expect("non-kernel point should avoid poles in y pullback");

    assert_eq!(
        Some(&x_value),
        crate::elliptic_curves::AffinePoint::x_coordinate(&image)
    );
    assert_eq!(
        Some(&y_value),
        crate::elliptic_curves::AffinePoint::y_coordinate(&image)
    );
}

#[test]
fn function_field_map_from_verschiebung_rejects_non_characteristic_scalar() {
    let curve = curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let candidate_v = crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
        frobenius.codomain().clone(),
        frobenius.domain().clone(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .x(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .y(),
    )
    .expect("identity candidate on the twist should validate");
    let verschiebung =
        VerschiebungIsogeny::new(frobenius.clone(), candidate_v).expect("candidate should build");
    let expected_left = frobenius
        .as_function_field_map()
        .compose(verschiebung.as_function_field_map())
        .expect("left composition should build");
    let expected_right = verschiebung
        .as_function_field_map()
        .compose(&frobenius.as_function_field_map())
        .expect("right composition should build");
    let certificate = VerschiebungCertificate::new(verschiebung, expected_left, expected_right)
        .expect("certificate should build");

    let scalar =
        ScalarMultiplicationIsogeny::new(curve, 2).expect("scalar multiplication should build");

    assert_eq!(
        scalar.as_function_field_map_from_verschiebung(&certificate),
        Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch))
    );
}

fn curve_and_point() -> impl Strategy<Value = (Curve, <Curve as CurveModel>::Point)> {
    let curve = curve();
    let points = curve.points();
    let len = points.len();

    (0usize..len).prop_map(move |index| (curve.clone(), points[index].clone()))
}

fn evaluate_short_weierstrass_function_at_point<F: Field>(
    function: &ShortWeierstrassFunction<F>,
    point: &crate::elliptic_curves::AffinePoint<F>,
) -> Option<F::Elem> {
    let x = crate::elliptic_curves::AffinePoint::x_coordinate(point)?;
    let y = crate::elliptic_curves::AffinePoint::y_coordinate(point)?;
    let a_value = evaluate_rational_function_at_x(function.a_part(), x)?;
    let b_value = evaluate_rational_function_at_x(function.b_part(), x)?;

    Some(F::add(&a_value, &F::mul(y, &b_value)))
}

fn evaluate_rational_function_at_x<F: Field>(
    function: &crate::fields::RationalFunction<F>,
    x: &F::Elem,
) -> Option<F::Elem> {
    let numerator = evaluate_dense(function.numerator(), x).ok()?;
    let denominator = evaluate_dense(function.denominator(), x).ok()?;

    if F::is_zero(&denominator) {
        None
    } else {
        F::div(&numerator, &denominator).ok()
    }
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

#[test]
fn characteristic_divisible_scalar_reports_mixed_kernel_description() {
    let curve = curve();
    let isogeny = ScalarMultiplicationIsogeny::new(curve, 41).expect("scalar isogeny should build");

    let description = isogeny.kernel_description();
    assert_eq!(description.reduced_degree(), Some(1));
    assert_eq!(description.infinitesimal_degree(), Some(41 * 41));
    assert_eq!(description.degree(), Some(41 * 41));
    assert_eq!(
        description.rational_points(),
        Some([isogeny.domain().identity()].as_slice())
    );
    assert!(matches!(description, KernelDescription::Mixed(_)));
}

#[test]
fn scalar_characteristic_factorization_splits_off_the_prime_to_p_part() {
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve(), 82).expect("scalar isogeny should build");
    let factorization = isogeny.scalar_characteristic_factorization();

    assert_eq!(factorization.p_power_exponent(), 1);
    assert_eq!(factorization.separable_part(), 2);
    assert_eq!(factorization.separable_degree(), 4);
    assert_eq!(factorization.infinitesimal_degree(), 41 * 41);
}

#[test]
fn visible_reduced_kernel_points_for_characteristic_divisible_scalar_come_from_the_prime_to_p_part()
{
    let curve = curve();
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 82).expect("scalar isogeny should build");

    let mut expected = vec![curve.identity()];
    expected.extend(curve.points_of_order(2));

    assert_eq!(
        isogeny
            .visible_reduced_kernel_points()
            .expect("visible reduced points should enumerate"),
        expected
    );
}

#[test]
fn mixed_kernel_description_for_p_times_m_uses_the_visible_m_torsion_and_residual_p_power_degree() {
    let curve = curve();
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 82).expect("scalar isogeny should build");
    let description = isogeny.kernel_description();

    let mut expected = vec![curve.identity()];
    expected.extend(curve.points_of_order(2));

    assert_eq!(description.reduced_degree(), Some(4));
    assert_eq!(description.infinitesimal_degree(), Some(41 * 41));
    assert_eq!(description.degree(), Some(4 * 41 * 41));
    assert_eq!(description.rational_points(), Some(expected.as_slice()));
    match description {
        KernelDescription::Mixed(mixed) => {
            assert_eq!(
                mixed.label(),
                Some("kernel contribution from [n] = [p^1] o [2]")
            );
        }
        other => panic!("expected mixed kernel description, got {other:?}"),
    }
}
