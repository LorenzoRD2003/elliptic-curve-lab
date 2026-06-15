use proptest::prelude::*;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusCharacteristicPolynomial,
        characteristic_equation::FrobeniusCharacteristicEquationCurveModel,
    },
    traits::{AffineCurveModel, EnumerableCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::{
    Fp,
    traits::{Field, FiniteField},
};
use crate::proptest_support::{config::CurveStrategyConfig, elliptic_curves::arb_curve_and_point};

type F41 = Fp<41>;
type F43 = Fp<43>;

#[test]
fn characteristic_equation_holds_on_a_prime_field_rational_point() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0,1) should lie on the curve");
    let characteristic_polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();

    let check = curve
        .verify_frobenius_characteristic_equation_at_point(&point, &characteristic_polynomial)
        .expect("characteristic equation should evaluate");
    assert!(check.holds());
}

#[test]
fn exhaustive_characteristic_equation_report_holds_on_a_small_prime_field_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let report = curve
        .verify_frobenius_characteristic_equation_exhaustive()
        .expect("exhaustive characteristic-equation check should evaluate");

    assert_eq!(report.checked_points(), curve.order());
    assert!(report.all_hold());
}

#[test]
fn characteristic_equation_rejects_incompatible_frobenius_base_field() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0,1) should lie on the curve");
    let incompatible_polynomial = FrobeniusCharacteristicPolynomial::new(
        crate::fields::finite_field_descriptor::FiniteFieldDescriptor::new(
            F41::characteristic(),
            F41::extension_degree(),
        )
        .expect("F41 descriptor should be valid"),
        0,
    );

    assert_eq!(
        curve.verify_frobenius_characteristic_equation_at_point(&point, &incompatible_polynomial),
        Err(CurveError::IncompatibleFrobeniusBaseField {
            curve_characteristic: F43::characteristic(),
            curve_extension_degree: F43::extension_degree().get(),
            polynomial_characteristic: F41::characteristic(),
            polynomial_extension_degree: F41::extension_degree().get(),
        })
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[test]
    fn property_characteristic_equation_holds_for_sampled_f43_rational_points(
        (curve, point) in arb_curve_and_point::<43>(CurveStrategyConfig::default()),
    ) {
        let characteristic_polynomial = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves")
            .characteristic_polynomial();

        let check = curve.verify_frobenius_characteristic_equation_at_point(&point, &characteristic_polynomial)
            .expect("characteristic equation should evaluate on enumerated rational points");

        prop_assert!(check.holds());
    }
}
