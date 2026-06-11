use crate::elliptic_curves::{
    CurveModel, EnumerableCurveModel, ShortWeierstrassCurve, absolute_frobenius_power_point,
    frobenius_twist_power, relative_frobenius_point,
};
use crate::fields::{Field, Fp, RationalFunction};
use crate::isogenies::{
    AbsoluteFrobeniusIsogeny, DegreeFactorizedIsogeny, FrobeniusLikeIsogeny, Isogeny,
    IsogenySeparabilityKind, RelativeFrobeniusIsogeny, VerschiebungCertificate, VerschiebungError,
    VerschiebungIsogeny,
};
use crate::polynomials::DensePolynomial;

type F17 = Fp<17>;

crate::fields::define_fp_quadratic_extension!(
    spec: F17Sqrt3FrobeniusIsoSpec,
    field: F17Sqrt3FrobeniusIso,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3)) for Frobenius isogeny tests",
);

fn extension_curve() -> ShortWeierstrassCurve<F17Sqrt3FrobeniusIso> {
    let generator = F17Sqrt3FrobeniusIso::element(vec![F17::zero(), F17::one()]);
    ShortWeierstrassCurve::<F17Sqrt3FrobeniusIso>::new(generator, F17Sqrt3FrobeniusIso::one())
        .expect("curve should be nonsingular")
}

fn prime_curve() -> ShortWeierstrassCurve<F17> {
    ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
        .expect("curve should be nonsingular")
}

fn extension_rhs(
    curve: &ShortWeierstrassCurve<F17Sqrt3FrobeniusIso>,
) -> RationalFunction<F17Sqrt3FrobeniusIso> {
    RationalFunction::<F17Sqrt3FrobeniusIso>::from_polynomial(
        DensePolynomial::<F17Sqrt3FrobeniusIso>::new(vec![
            curve.b().clone(),
            curve.a().clone(),
            F17Sqrt3FrobeniusIso::zero(),
            F17Sqrt3FrobeniusIso::one(),
        ]),
    )
}

#[test]
fn absolute_frobenius_isogeny_uses_the_p_power_twist_and_degree_factorization() {
    let curve = extension_curve();
    let isogeny =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");

    assert_eq!(
        isogeny.codomain(),
        &frobenius_twist_power(&curve, 1).expect("twist should build")
    );
    assert_eq!(isogeny.separable_degree(), 1);
    assert_eq!(isogeny.inseparable_degree(), 17);
    assert!(isogeny.is_purely_inseparable());
    assert!(isogeny.kernel_points().is_empty());
}

#[test]
fn absolute_frobenius_point_evaluation_matches_existing_helper() {
    let curve = extension_curve();
    let isogeny =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| curve.is_on_curve_nonzero(point))
        .expect("curve should have a finite rational point");

    assert_eq!(
        isogeny.evaluate(&point).expect("evaluation should work"),
        absolute_frobenius_power_point(&curve, &point, 1).expect("existing helper should work")
    );
}

#[test]
fn absolute_frobenius_pullbacks_match_x_to_the_p_and_y_to_the_p() {
    let curve = extension_curve();
    let isogeny =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let rhs = extension_rhs(&curve);

    assert_eq!(
        isogeny.x_pullback().a_part(),
        &RationalFunction::<F17Sqrt3FrobeniusIso>::from_polynomial(DensePolynomial::<
            F17Sqrt3FrobeniusIso,
        >::new({
            let mut coefficients = vec![F17Sqrt3FrobeniusIso::zero(); 18];
            coefficients[17] = F17Sqrt3FrobeniusIso::one();
            coefficients
        }))
    );
    assert!(isogeny.x_pullback().b_part().is_zero());

    let mut rhs_to_the_eighth =
        RationalFunction::<F17Sqrt3FrobeniusIso>::constant(F17Sqrt3FrobeniusIso::one());
    for _ in 0..8 {
        rhs_to_the_eighth = rhs_to_the_eighth.mul(&rhs);
    }
    assert!(isogeny.y_pullback().a_part().is_zero());
    assert_eq!(isogeny.y_pullback().b_part(), &rhs_to_the_eighth);
}

#[test]
fn relative_frobenius_isogeny_is_pointwise_identity_on_the_represented_field() {
    let curve = extension_curve();
    let isogeny =
        RelativeFrobeniusIsogeny::new(curve.clone()).expect("relative Frobenius should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| curve.is_on_curve_nonzero(point))
        .expect("curve should have a finite rational point");

    assert_eq!(isogeny.domain(), isogeny.codomain());
    assert_eq!(isogeny.separable_degree(), 1);
    assert_eq!(isogeny.inseparable_degree(), 17_u128.pow(2));
    assert_eq!(
        isogeny.evaluate(&point).expect("evaluation should work"),
        relative_frobenius_point(&curve, &point).expect("existing helper should work")
    );
}

#[test]
fn frobenius_differential_reports_are_classified_as_purely_inseparable() {
    let absolute =
        AbsoluteFrobeniusIsogeny::new(extension_curve()).expect("absolute Frobenius should build");
    let relative =
        RelativeFrobeniusIsogeny::new(prime_curve()).expect("relative Frobenius should build");

    let absolute_report = absolute
        .differential_pullback_report()
        .expect("absolute report should build");
    let relative_report = relative
        .differential_pullback_report()
        .expect("relative report should build");

    assert!(absolute_report.dx_pullback().is_zero());
    assert!(relative_report.dx_pullback().is_zero());
    assert_eq!(
        absolute_report.separability_kind(),
        IsogenySeparabilityKind::PurelyInseparable
    );
    assert_eq!(
        relative_report.separability_kind(),
        IsogenySeparabilityKind::PurelyInseparable
    );
    assert!(
        absolute_report
            .invariant_differential_multiplier()
            .is_zero()
    );
    assert!(
        relative_report
            .invariant_differential_multiplier()
            .is_zero()
    );
}

#[test]
fn frobenius_like_trait_supplies_shared_pullback_methods_without_an_enum_wrapper() {
    let absolute =
        AbsoluteFrobeniusIsogeny::new(prime_curve()).expect("absolute Frobenius should build");
    let relative =
        RelativeFrobeniusIsogeny::new(prime_curve()).expect("relative Frobenius should build");

    assert_eq!(absolute.separable_degree(), 1);
    assert_eq!(absolute.inseparable_degree(), 17);
    assert_eq!(relative.separable_degree(), 1);
    assert_eq!(relative.inseparable_degree(), 17);
    assert!(
        absolute
            .as_function_field_map()
            .differential_pullback_report()
            .is_ok()
    );
    assert!(
        relative
            .as_function_field_map()
            .differential_pullback_report()
            .is_ok()
    );
}

#[test]
fn verschiebung_constructor_rejects_mismatched_pullback_direction() {
    let curve = extension_curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let field = crate::elliptic_curves::ShortWeierstrassFunctionField::<F17Sqrt3FrobeniusIso>::new(
        curve.clone(),
    );
    let wrong_direction = crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve,
        field.x(),
        field.y(),
    )
    .expect("identity map should validate");

    assert!(matches!(
        VerschiebungIsogeny::new(frobenius, wrong_direction),
        Err(crate::isogenies::IsogenyError::Verschiebung(
            VerschiebungError::DomainCodomainMismatch
        ))
    ));
}

#[test]
fn verschiebung_verification_uses_pullback_composition_relations() {
    let curve = prime_curve();
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

    assert_eq!(verschiebung.degree(), 17);
    assert_eq!(verschiebung.domain_curve(), frobenius.codomain());
    assert_eq!(verschiebung.codomain_curve(), frobenius.domain());
    assert_eq!(
        verschiebung.verify_v_after_f_equals_p(&expected_left),
        Ok(())
    );
    assert_eq!(
        verschiebung.verify_f_after_v_equals_p(&expected_right),
        Ok(())
    );
    assert_eq!(
        verschiebung.verify_duality_relations(&expected_left, &expected_right),
        Ok(())
    );
}

#[test]
fn verschiebung_verification_rejects_wrong_expected_pullback() {
    let curve = prime_curve();
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
        VerschiebungIsogeny::new(frobenius, candidate_v).expect("candidate should build");
    let wrong_expected = crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve.clone(),
        crate::elliptic_curves::ShortWeierstrassFunctionField::<F17>::new(curve.clone()).x(),
        crate::elliptic_curves::ShortWeierstrassFunctionField::<F17>::new(curve).y(),
    )
    .expect("identity map should validate");

    assert_eq!(
        verschiebung.verify_v_after_f_equals_p(&wrong_expected),
        Err(crate::isogenies::IsogenyError::Verschiebung(
            VerschiebungError::LeftDualityViolation
        ))
    );
}

#[test]
fn verschiebung_certificate_packages_the_expected_maps_and_verifies_without_arguments() {
    let curve = prime_curve();
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
        VerschiebungCertificate::new(verschiebung, expected_left.clone(), expected_right.clone())
            .expect("certificate should build");

    assert_eq!(certificate.frobenius().domain(), frobenius.domain());
    assert_eq!(certificate.multiplication_by_p_on_e(), &expected_left);
    assert_eq!(
        certificate.multiplication_by_p_on_frobenius_twist(),
        &expected_right
    );
    assert_eq!(certificate.verify_v_after_f_equals_p(), Ok(()));
    assert_eq!(certificate.verify_f_after_v_equals_p(), Ok(()));
    assert_eq!(certificate.verify_duality_relations(), Ok(()));
}

#[test]
fn verschiebung_candidate_can_upgrade_directly_to_a_certificate() {
    let curve = prime_curve();
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

    let certificate = verschiebung
        .certify(expected_left, expected_right)
        .expect("candidate should certify");

    assert_eq!(certificate.verify_duality_relations(), Ok(()));
}
