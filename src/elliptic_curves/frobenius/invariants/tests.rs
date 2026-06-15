use core::num::NonZeroU32;
use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{FrobeniusCurveType, FrobeniusTrace},
    traits::FrobeniusTraceCurveModel,
};
use crate::fields::{Fp, finite_field_descriptor::FiniteFieldDescriptor, traits::Field};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};

type F43 = Fp<43>;

fn nz(n: u32) -> NonZeroU32 {
    NonZeroU32::new(n).expect("test degrees are positive")
}

#[test]
fn frobenius_trace_from_order_and_order_from_trace_roundtrip() {
    let base_field = FiniteFieldDescriptor::new(43, nz(1)).expect("F43 descriptor should be valid");
    let report = FrobeniusTrace::from_order(base_field.clone(), 48)
        .expect("small Frobenius trace package should build");

    assert_eq!(report.trace(), -4);
    assert_eq!(
        FrobeniusTrace::curve_order_from_trace(base_field, report.trace()),
        Ok(report.curve_order())
    );
}

#[test]
fn characteristic_polynomial_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let polynomial = trace.characteristic_polynomial();

    assert_eq!(polynomial.base_field(), trace.base_field());
    assert_eq!(polynomial.trace(), trace.trace());
    assert_eq!(polynomial.field_order(), trace.field_order());
}

#[test]
fn frobenius_discriminant_matches_the_quadratic_formula() {
    let base_field = FiniteFieldDescriptor::new(43, nz(1)).expect("F43 descriptor should be valid");
    let trace = FrobeniusTrace::from_order(base_field, 41).expect("t = 3 should be valid over F43");
    let discriminant = trace.discriminant();

    assert_eq!(
        discriminant.quadratic_discriminant().value(),
        &num_bigint::BigInt::from(-163)
    );
    assert!(discriminant.is_negative());
    assert!(discriminant.is_fundamental());
}

#[test]
fn local_zeta_function_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let zeta = trace.local_zeta_function();

    assert_eq!(
        zeta,
        trace.characteristic_polynomial().local_zeta_function()
    );
}

#[test]
fn curve_type_report_classifies_ordinary_and_supersingular_examples() {
    let ordinary = FrobeniusTrace::from_order(
        FiniteFieldDescriptor::new(43, nz(1)).expect("F43 descriptor should be valid"),
        48,
    )
    .expect("t = -4 should be valid");
    let supersingular = FrobeniusTrace::from_order(
        FiniteFieldDescriptor::new(43, nz(1)).expect("F43 descriptor should be valid"),
        44,
    )
    .expect("t = 0 should be valid");

    assert_eq!(ordinary.curve_type(), FrobeniusCurveType::Ordinary);
    assert_eq!(
        supersingular.curve_type(),
        FrobeniusCurveType::Supersingular
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]

    #[test]
    fn property_characteristic_polynomial_and_zeta_are_consistent_with_trace(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let polynomial = trace.characteristic_polynomial();
        let zeta_from_polynomial = polynomial.local_zeta_function();
        let zeta_from_trace = trace.local_zeta_function();

        prop_assert_eq!(polynomial.base_field(), trace.base_field());
        prop_assert_eq!(polynomial.trace(), trace.trace());
        prop_assert_eq!(polynomial.field_order(), trace.field_order());
        prop_assert_eq!(polynomial.evaluate_at_integer(1), i128::from(trace.curve_order()));
        prop_assert_eq!(&zeta_from_polynomial, &zeta_from_trace);
    }
}
