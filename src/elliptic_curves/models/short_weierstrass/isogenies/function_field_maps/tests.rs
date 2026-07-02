use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::function_fields::ShortWeierstrassFunction,
    short_weierstrass::isogenies::function_field_maps::{
        IsogenySeparabilityKind, ShortWeierstrassFunctionFieldMap,
    },
};
use crate::fields::traits::*;
use crate::fields::{Q, rational_function_field::RationalFunction};
use crate::isogenies::error::{IsogenyError, IsogenyMapError};
use crate::polynomials::DensePolynomial;

type F17 = crate::fields::Fp17;

fn f17_curve() -> ShortWeierstrassCurve<F17> {
    ShortWeierstrassCurve::new(F17::from_i64(2), F17::from_i64(3)).expect("valid curve")
}

fn alternate_f17_curve() -> ShortWeierstrassCurve<F17> {
    ShortWeierstrassCurve::new(F17::from_i64(0), F17::from_i64(1)).expect("valid curve")
}

fn q_curve() -> ShortWeierstrassCurve<Q> {
    ShortWeierstrassCurve::new(Q::from_i64(-1), Q::zero()).expect("valid curve")
}

fn f17_dense(coefficients: &[i64]) -> DensePolynomial<F17> {
    DensePolynomial::new(coefficients.iter().copied().map(F17::from_i64).collect())
}

fn q_dense(coefficients: &[i64]) -> DensePolynomial<Q> {
    DensePolynomial::new(coefficients.iter().copied().map(Q::from_i64).collect())
}

#[test]
fn identity_pullback_is_valid_and_exposes_curve_accessors() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let map =
        ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve.clone(), field.x(), field.y())
            .expect("identity pullback should validate");

    assert_eq!(map.domain_curve(), &curve);
    assert_eq!(map.codomain_curve(), &curve);
    assert_eq!(map.x_pullback(), &field.x());
    assert_eq!(map.y_pullback(), &field.y());
    assert_eq!(map.domain_function_field().x(), field.x());
    assert_eq!(map.codomain_function_field().y(), field.y());
}

#[test]
fn constructor_rejects_pullbacks_that_do_not_live_on_the_declared_domain_curve() {
    let domain = f17_curve();
    let wrong_curve = alternate_f17_curve();
    let wrong_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(wrong_curve.clone());

    let result = ShortWeierstrassFunctionFieldMap::new(
        domain.clone(),
        domain,
        wrong_field.x(),
        wrong_field.y(),
    );

    assert_eq!(
        result,
        Err(IsogenyError::Map(
            IsogenyMapError::FunctionFieldMapPullbackCurveMismatch
        ))
    );
}

#[test]
fn constructor_rejects_pullbacks_that_do_not_satisfy_the_codomain_equation() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());

    let result = ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve, field.x(), field.x());

    assert_eq!(
        result,
        Err(IsogenyError::Map(
            IsogenyMapError::FunctionFieldMapCodomainEquationViolation
        ))
    );
}

#[test]
fn pullback_of_a_rational_function_substitutes_x_pullback() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let map =
        ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve.clone(), field.x(), field.y())
            .expect("identity pullback should validate");
    let rational = RationalFunction::<F17>::new(f17_dense(&[1, 0, 1]), f17_dense(&[1, 1]))
        .expect("denominator should be non-zero");

    let pulled = map
        .pullback_rational_function(&rational)
        .expect("identity pullback should preserve rational functions");

    assert_eq!(
        pulled,
        ShortWeierstrassFunction::from_rational_function(curve, rational)
    );
}

#[test]
fn pullback_of_a_function_uses_the_stored_y_pullback() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let negation = ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve.clone(),
        field.x(),
        field.y().neg(),
    )
    .expect("negation pullback should validate");
    let function = ShortWeierstrassFunction::new(
        curve.clone(),
        RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1]))
            .expect("denominator should be non-zero"),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[2, 1])),
    );
    let expected =
        ShortWeierstrassFunction::new(curve, function.a_part().clone(), function.b_part().neg());

    assert_eq!(
        negation
            .pullback_function(&function)
            .expect("negation pullback should evaluate"),
        expected
    );
}

#[test]
fn pullback_function_rejects_functions_from_the_wrong_codomain_curve() {
    let domain = f17_curve();
    let codomain = alternate_f17_curve();
    let domain_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(domain.clone());
    let map = ShortWeierstrassFunctionFieldMap::new(
        domain.clone(),
        codomain,
        ShortWeierstrassFunction::<F17>::from_rational_function(
            domain.clone(),
            RationalFunction::<F17>::constant(F17::zero()),
        ),
        ShortWeierstrassFunction::<F17>::from_rational_function(
            domain.clone(),
            RationalFunction::<F17>::constant(F17::one()),
        ),
    )
    .expect("constant codomain point should satisfy the codomain equation");

    assert_eq!(
        map.pullback_function(&domain_field.x()),
        Err(IsogenyError::Map(
            IsogenyMapError::FunctionFieldMapSourceCurveMismatch
        ))
    );
}

#[test]
fn pullback_rational_function_reports_when_the_denominator_maps_to_zero() {
    let curve = q_curve();
    let constant_point_map = ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve.clone(),
        ShortWeierstrassFunction::<Q>::from_rational_function(
            curve.clone(),
            RationalFunction::<Q>::constant(Q::zero()),
        ),
        ShortWeierstrassFunction::<Q>::from_rational_function(
            curve.clone(),
            RationalFunction::<Q>::constant(Q::zero()),
        ),
    )
    .expect("a constant rational point still satisfies the codomain equation");
    let rational = RationalFunction::<Q>::new(q_dense(&[1]), q_dense(&[0, 1]))
        .expect("denominator should be non-zero");

    assert_eq!(
        constant_point_map.pullback_rational_function(&rational),
        Err(IsogenyError::Map(
            IsogenyMapError::FunctionFieldMapDenominatorMapsToZero
        ))
    );
}

#[test]
fn composition_is_contravariant_on_pullbacks() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let negation = ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve.clone(),
        field.x(),
        field.y().neg(),
    )
    .expect("negation pullback should validate");
    let identity =
        ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve.clone(), field.x(), field.y())
            .expect("identity pullback should validate");

    assert_eq!(
        negation
            .compose(&negation)
            .expect("negation composed with negation should validate"),
        identity
    );
}

#[test]
fn composition_rejects_mismatched_middle_curves() {
    let first_curve = f17_curve();
    let middle_curve = alternate_f17_curve();
    let last_curve =
        ShortWeierstrassCurve::new(F17::from_i64(6), F17::from_i64(4)).expect("valid curve");
    let last_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(last_curve.clone());
    let first = ShortWeierstrassFunctionFieldMap::new(
        first_curve.clone(),
        middle_curve.clone(),
        ShortWeierstrassFunction::<F17>::from_rational_function(
            first_curve.clone(),
            RationalFunction::<F17>::constant(F17::zero()),
        ),
        ShortWeierstrassFunction::<F17>::from_rational_function(
            first_curve.clone(),
            RationalFunction::<F17>::constant(F17::one()),
        ),
    )
    .expect("constant codomain point should satisfy the middle curve equation");
    let second = ShortWeierstrassFunctionFieldMap::new(
        last_curve.clone(),
        last_curve.clone(),
        last_field.x(),
        last_field.y(),
    )
    .expect("second identity pullback should validate");

    assert_eq!(
        first.compose(&second),
        Err(IsogenyError::Map(
            IsogenyMapError::CompositionDomainCodomainMismatch
        ))
    );
}

#[test]
fn identity_differential_pullback_report_has_unit_multiplier_and_separable_kind() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let map =
        ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve.clone(), field.x(), field.y())
            .expect("identity pullback should validate");
    let report = map
        .differential_pullback_report()
        .expect("identity report should build");

    assert!(report.dx_pullback().is_one());
    assert_eq!(
        report
            .rational_multiplier()
            .expect("identity multiplier should lie in F(x)"),
        &RationalFunction::<F17>::constant(F17::one())
    );
    assert_eq!(
        report.separability_kind(),
        IsogenySeparabilityKind::Separable
    );
    assert!(report.is_certified_separable());
}

#[test]
fn negation_differential_pullback_report_has_minus_one_multiplier_and_separable_kind() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let map = ShortWeierstrassFunctionFieldMap::new(
        curve.clone(),
        curve.clone(),
        field.x(),
        field.y().neg(),
    )
    .expect("negation pullback should validate");
    let report = map
        .differential_pullback_report()
        .expect("negation report should build");

    assert_eq!(
        report
            .rational_multiplier()
            .expect("negation multiplier should lie in F(x)"),
        &RationalFunction::<F17>::constant(F17::from_i64(-1))
    );
    assert_eq!(
        report.separability_kind(),
        IsogenySeparabilityKind::Separable
    );
}

#[test]
fn constant_maps_are_currently_classified_as_constant_or_invalid() {
    let curve = alternate_f17_curve();
    let domain = f17_curve();
    let map = ShortWeierstrassFunctionFieldMap::new(
        domain.clone(),
        curve.clone(),
        ShortWeierstrassFunction::<F17>::from_rational_function(
            domain.clone(),
            RationalFunction::<F17>::constant(F17::zero()),
        ),
        ShortWeierstrassFunction::<F17>::from_rational_function(
            domain,
            RationalFunction::<F17>::constant(F17::one()),
        ),
    )
    .expect("constant codomain point should satisfy the codomain equation");
    let report = map
        .differential_pullback_report()
        .expect("constant report should build when Y is non-zero");

    assert!(report.invariant_differential_multiplier().is_zero());
    assert_eq!(
        report.separability_kind(),
        IsogenySeparabilityKind::ConstantOrInvalid
    );
    assert!(!report.is_certified_separable());
}

#[test]
fn frobenius_style_nonconstant_zero_dx_map_is_classified_as_purely_inseparable() {
    let curve = f17_curve();
    let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));

    let x_pullback = ShortWeierstrassFunction::<F17>::from_rational_function(
        curve.clone(),
        RationalFunction::<F17>::from_polynomial(DensePolynomial::<F17>::new({
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[17] = F17::one();
            coefficients
        })),
    );

    let mut rhs_to_the_eighth = RationalFunction::<F17>::constant(F17::one());
    for _ in 0..8 {
        rhs_to_the_eighth = rhs_to_the_eighth.mul(&rhs);
    }
    let y_pullback = ShortWeierstrassFunction::<F17>::new(
        curve.clone(),
        RationalFunction::<F17>::constant(F17::zero()),
        rhs_to_the_eighth,
    );

    let map = ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve, x_pullback, y_pullback)
        .expect("frobenius-style pullback should validate");
    let report = map
        .differential_pullback_report()
        .expect("report should build");

    assert!(report.dx_pullback().is_zero());
    assert!(report.invariant_differential_multiplier().is_zero());
    assert_eq!(
        report.separability_kind(),
        IsogenySeparabilityKind::PurelyInseparable
    );
}
