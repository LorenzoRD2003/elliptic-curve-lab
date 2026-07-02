use crate::fields::traits::*;
use proptest::prelude::*;

use crate::fields::{Q, traits::PthRootExtraction};
use crate::polynomials::{
    DensePolynomial, PolynomialError, SparsePolynomial, sparse::SparsePolynomialTerm,
    traits::UnivariatePolynomial,
};
use crate::proptest_support::{
    config::PolynomialStrategyConfig, polynomials::arb_dense_polynomial,
};

type F17 = crate::fields::Fp17;
type F17Samples = Vec<(
    <F17 as crate::fields::traits::Field>::Elem,
    <F17 as crate::fields::traits::Field>::Elem,
)>;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F17Sqrt3DensePthRootSpec,
    field: F17Sqrt3DensePthRoot,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3)) for dense polynomial p-th-root tests",
);

fn f17_coefficients(values: &[u64]) -> Vec<<F17 as crate::fields::traits::Field>::Elem> {
    values.iter().copied().map(F17::from_i64).collect()
}

fn q_coefficients(values: &[(i64, i64)]) -> Vec<<Q as crate::fields::traits::Field>::Elem> {
    values
        .iter()
        .map(|&(numerator, denominator)| {
            let numerator = Q::from_i64(numerator);
            let denominator = Q::from_i64(denominator);
            Q::div(&numerator, &denominator).expect("denominator should be non-zero")
        })
        .collect()
}

fn q(numerator: i64, denominator: i64) -> <Q as crate::fields::traits::Field>::Elem {
    let numerator = Q::from_i64(numerator);
    let denominator = Q::from_i64(denominator);
    Q::div(&numerator, &denominator).expect("denominator should be non-zero")
}

fn assert_dense_eq<F: Field>(actual: &DensePolynomial<F>, expected: &DensePolynomial<F>) {
    assert_eq!(actual.coefficients().len(), expected.coefficients().len());

    for (actual, expected) in actual.coefficients().iter().zip(expected.coefficients()) {
        assert!(F::eq(actual, expected));
    }
}

#[test]
fn evaluate_dense_uses_horner_rule_over_f17() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 2]));
    let value = polynomial
        .evaluate(&F17::from_i64(4))
        .expect("evaluation should work");

    assert!(F17::eq(&value, &F17::from_i64(4)));
}

#[test]
fn evaluate_dense_zero_polynomial_returns_zero() {
    let polynomial = DensePolynomial::<F17>::new(Vec::new());
    let value = polynomial
        .evaluate(&F17::from_i64(9))
        .expect("evaluation should work");

    assert!(F17::eq(&value, &F17::zero()));
}

#[test]
fn lagrange_interpolate_returns_zero_for_empty_input() {
    let polynomial =
        DensePolynomial::<F17>::lagrange_interpolate(&[]).expect("empty interpolation should work");

    assert!(polynomial.is_zero());
    assert_eq!(polynomial.coefficients(), &[]);
}

#[test]
fn lagrange_interpolate_returns_constant_for_single_sample() {
    let polynomial =
        DensePolynomial::<F17>::lagrange_interpolate(&[(F17::from_i64(9), F17::from_i64(4))])
            .expect("single sample should interpolate");

    assert_dense_eq(
        &polynomial,
        &DensePolynomial::<F17>::new(vec![F17::from_i64(4)]),
    );
}

#[test]
fn lagrange_interpolate_reconstructs_linear_polynomial_over_f17() {
    let samples = [
        (F17::from_i64(0), F17::from_i64(3)),
        (F17::from_i64(1), F17::from_i64(8)),
    ];

    let polynomial =
        DensePolynomial::<F17>::lagrange_interpolate(&samples).expect("interpolation should work");

    assert_dense_eq(
        &polynomial,
        &DensePolynomial::<F17>::new(vec![F17::from_i64(3), F17::from_i64(5)]),
    );
}

#[test]
fn lagrange_interpolate_reconstructs_quadratic_polynomial_over_f17() {
    let samples = [
        (F17::from_i64(0), F17::from_i64(3)),
        (F17::from_i64(1), F17::from_i64(10)),
        (F17::from_i64(2), F17::from_i64(4)),
    ];

    let polynomial =
        DensePolynomial::<F17>::lagrange_interpolate(&samples).expect("interpolation should work");

    assert_dense_eq(
        &polynomial,
        &DensePolynomial::<F17>::new(vec![F17::from_i64(3), F17::from_i64(5), F17::from_i64(2)]),
    );
}

#[test]
fn lagrange_interpolate_matches_all_input_samples_over_q() {
    let samples = [(q(0, 1), q(1, 2)), (q(1, 1), q(7, 6)), (q(2, 1), q(17, 6))];

    let polynomial =
        DensePolynomial::<Q>::lagrange_interpolate(&samples).expect("interpolation should work");

    assert_dense_eq(
        &polynomial,
        &DensePolynomial::<Q>::new(vec![q(1, 2), q(1, 6), q(1, 2)]),
    );

    for (x, y) in &samples {
        let value = polynomial.evaluate(x).expect("evaluation should work");
        assert!(Q::eq(&value, y));
    }
}

#[test]
fn lagrange_interpolate_rejects_duplicate_x_coordinates() {
    let samples = [
        (F17::from_i64(3), F17::from_i64(1)),
        (F17::from_i64(3), F17::from_i64(9)),
    ];

    let error = DensePolynomial::<F17>::lagrange_interpolate(&samples)
        .expect_err("duplicate x values should fail");

    assert_eq!(error, PolynomialError::DuplicateInterpolationAbscissa);
}

fn interpolation_case() -> impl Strategy<Value = (DensePolynomial<F17>, F17Samples)> {
    arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
        max_len: 4,
        ..PolynomialStrategyConfig::default()
    })
    .prop_flat_map(|polynomial| {
        let sample_count = polynomial.degree().map_or(1, |degree| degree + 1);
        crate::proptest_support::fields::arb_distinct_fp_elems::<crate::fields::Fp17>(sample_count)
            .prop_map(move |xs| {
                let samples = xs
                    .into_iter()
                    .map(|x| {
                        let y = polynomial.evaluate(&x).expect("evaluation should succeed");
                        (x, y)
                    })
                    .collect::<Vec<_>>();
                (polynomial.clone(), samples)
            })
    })
}

#[test]
fn dense_polynomial_preserves_storage_order_after_normalization() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[3, 15, 0, 7]));

    let coefficients = polynomial.coefficients();
    assert_eq!(coefficients.len(), 4);
    assert!(F17::eq(&coefficients[0], &F17::from_i64(3)));
    assert!(F17::eq(&coefficients[1], &F17::from_i64(15)));
    assert!(F17::eq(&coefficients[2], &F17::from_i64(0)));
    assert!(F17::eq(&coefficients[3], &F17::from_i64(7)));
    assert_eq!(polynomial.len(), 4);
    assert_eq!(polynomial.degree(), Some(3));
    assert!(F17::eq(
        polynomial
            .leading_coefficient()
            .expect("leading coefficient"),
        &F17::from_i64(7)
    ));
}

#[test]
fn dense_polynomial_allows_empty_storage_for_zero_representation() {
    let polynomial = DensePolynomial::<F17>::new(Vec::new());

    assert!(polynomial.is_zero());
    assert_eq!(polynomial.len(), 0);
    assert_eq!(polynomial.degree(), None);
    assert_eq!(polynomial.leading_coefficient(), None);
    assert_eq!(polynomial.constant_term(), None);
}

#[test]
fn dense_polynomial_trims_trailing_zero_coefficients() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[5, 0, 0]));

    assert_eq!(polynomial.coefficients().len(), 1);
    assert_eq!(polynomial.degree(), Some(0));
    assert!(F17::eq(
        polynomial
            .leading_coefficient()
            .expect("leading coefficient"),
        &F17::from_i64(5)
    ));
}

#[test]
fn dense_polynomial_normalizes_all_zero_storage_to_empty() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[0, 0, 0]));

    assert!(polynomial.is_zero());
    assert_eq!(polynomial.coefficients(), &[]);
    assert_eq!(polynomial.degree(), None);
    assert_eq!(polynomial.leading_coefficient(), None);
}

#[test]
fn dense_polynomial_single_coefficient_has_degree_zero() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[8]));

    assert_eq!(polynomial.degree(), Some(0));
    assert!(F17::eq(
        polynomial.constant_term().expect("constant term"),
        &F17::from_i64(8)
    ));
    assert!(F17::eq(
        polynomial
            .leading_coefficient()
            .expect("leading coefficient"),
        &F17::from_i64(8)
    ));
    assert!(!polynomial.is_zero());
}

#[test]
fn dense_polynomial_addition_is_coefficient_wise_over_f17() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[15, 14]));
    let sum = lhs.add(&rhs);

    let coefficients = sum.coefficients();
    assert_eq!(coefficients.len(), 3);
    assert!(F17::eq(&coefficients[0], &F17::from_i64(1)));
    assert!(F17::eq(&coefficients[1], &F17::from_i64(2)));
    assert!(F17::eq(&coefficients[2], &F17::from_i64(1)));
}

#[test]
fn dense_polynomial_pth_root_over_prime_field_recovers_expected_coefficients() {
    let mut coefficients = vec![F17::zero(); 35];
    coefficients[0] = F17::from_i64(4);
    coefficients[17] = F17::from_i64(9);
    coefficients[34] = F17::from_i64(3);
    let polynomial = DensePolynomial::<F17>::new(coefficients);

    let root = polynomial
        .pth_root()
        .expect("all non-zero term degrees are divisible by the characteristic");

    assert_eq!(
        root,
        DensePolynomial::<F17>::new(f17_coefficients(&[4, 9, 3]))
    );
    assert!(polynomial.has_pth_root());
}

#[test]
fn dense_polynomial_pth_root_rejects_non_multiple_degree_terms() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 0, 0, 3]));

    assert_eq!(polynomial.pth_root(), None);
    assert!(!polynomial.has_pth_root());
}

#[test]
fn dense_polynomial_pth_root_handles_the_zero_polynomial() {
    let polynomial = DensePolynomial::<F17>::new(Vec::new());

    assert_eq!(
        polynomial.pth_root(),
        Some(DensePolynomial::<F17>::new(Vec::new()))
    );
    assert!(polynomial.has_pth_root());
}

#[test]
fn dense_polynomial_pth_root_uses_extension_field_coefficient_roots() {
    let generator = F17Sqrt3DensePthRoot::element(vec![F17::zero(), F17::one()]);
    let expected_root = DensePolynomial::<F17Sqrt3DensePthRoot>::new(vec![
        generator.clone(),
        F17Sqrt3DensePthRoot::one(),
    ]);

    let mut coefficients = vec![F17Sqrt3DensePthRoot::zero(); 18];
    coefficients[0] = F17Sqrt3DensePthRoot::pow(
        &generator,
        &F17Sqrt3DensePthRoot::characteristic()
            .to_positive_biguint()
            .expect("finite fields should have positive characteristic"),
    );
    coefficients[17] = F17Sqrt3DensePthRoot::one();
    let polynomial = DensePolynomial::<F17Sqrt3DensePthRoot>::new(coefficients);

    assert!(polynomial.pth_root() == Some(expected_root));
    assert!(
        polynomial
            .constant_term()
            .expect("constant term should be present")
            .pth_root()
            == Some(generator)
    );
}

#[test]
fn dense_polynomial_multiplication_uses_naive_convolution_over_f17() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 4, 5]));
    let product = lhs.mul(&rhs);

    let coefficients = product.coefficients();
    assert_eq!(coefficients.len(), 4);
    assert!(F17::eq(&coefficients[0], &F17::from_i64(3)));
    assert!(F17::eq(&coefficients[1], &F17::from_i64(10)));
    assert!(F17::eq(&coefficients[2], &F17::from_i64(13)));
    assert!(F17::eq(&coefficients[3], &F17::from_i64(10)));
}

#[test]
fn dense_polynomial_addition_works_over_q_too() {
    let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3)]));
    let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 3), (-2, 3), (5, 4)]));
    let sum = lhs.add(&rhs);

    let coefficients = sum.coefficients();
    assert_eq!(coefficients.len(), 3);
    assert!(Q::eq(
        &coefficients[0],
        &Q::div(&Q::from_i64(5), &Q::from_i64(6)).unwrap()
    ));
    assert!(Q::eq(&coefficients[1], &Q::zero()));
    assert!(Q::eq(
        &coefficients[2],
        &Q::div(&Q::from_i64(5), &Q::from_i64(4)).unwrap()
    ));
}

#[test]
fn dense_polynomial_constant_constructor_is_canonical() {
    let polynomial = DensePolynomial::<F17>::constant(F17::from_i64(9));

    assert_eq!(polynomial.coefficients().len(), 1);
    assert!(F17::eq(
        polynomial.constant_term().expect("constant term"),
        &F17::from_i64(9)
    ));
}

#[test]
fn dense_polynomial_negation_and_subtraction_work_over_f17() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[15, 14]));

    let neg_rhs = rhs.neg();
    assert_eq!(
        neg_rhs,
        DensePolynomial::<F17>::new(f17_coefficients(&[2, 3]))
    );

    let difference = lhs.sub(&rhs);
    assert_eq!(
        difference,
        DensePolynomial::<F17>::new(f17_coefficients(&[5, 8, 1]))
    );
}

#[test]
fn dense_polynomial_scale_multiplies_every_coefficient() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
    let scaled = polynomial.scale(&F17::from_i64(4));

    assert_eq!(
        scaled,
        DensePolynomial::<F17>::new(f17_coefficients(&[12, 3, 4]))
    );
}

#[test]
fn dense_polynomial_derivative_drops_the_constant_term() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[4, 3, 5, 2]));

    assert_eq!(
        polynomial.derivative(),
        DensePolynomial::<F17>::new(f17_coefficients(&[3, 10, 6]))
    );
}

#[test]
fn dense_polynomial_derivative_of_constant_is_zero() {
    let polynomial = DensePolynomial::<F17>::constant(F17::from_i64(9));

    assert!(polynomial.derivative().is_zero());
    assert_eq!(polynomial.derivative().coefficients(), &[]);
}

#[test]
fn dense_polynomial_derivative_trims_zero_tail_after_characteristic_cancellation() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 0, 0, 0, 0, 17]));

    assert!(polynomial.derivative().is_zero());
    assert_eq!(polynomial.derivative().coefficients(), &[]);
}

#[test]
fn dense_polynomial_derivative_works_over_q_too() {
    let polynomial = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3), (3, 4)]));

    assert_eq!(
        polynomial.derivative(),
        DensePolynomial::<Q>::new(q_coefficients(&[(2, 3), (3, 2)]))
    );
}

#[test]
fn dense_polynomial_manual_partial_eq_uses_field_equality() {
    let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(2, 4), (3, 6)]));
    let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (1, 2)]));
    let different = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3)]));

    assert_eq!(lhs, rhs);
    assert_ne!(lhs, different);
}

#[test]
fn dense_polynomial_monic_helpers_work_over_f17() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[6, 3]));
    let monic = polynomial
        .make_monic()
        .expect("non-zero polynomial should normalize");

    assert!(!polynomial.is_monic());
    assert!(monic.is_monic());
    assert_eq!(
        monic,
        DensePolynomial::<F17>::new(f17_coefficients(&[2, 1]))
    );
}

#[test]
fn dense_polynomial_zero_has_no_monic_normalization() {
    let polynomial = DensePolynomial::<F17>::new(Vec::new());
    let error = polynomial
        .make_monic()
        .expect_err("zero polynomial should not be monic-normalizable");

    assert_eq!(
        error,
        PolynomialError::ZeroPolynomialHasNoMonicNormalization
    );
}

#[test]
fn dense_polynomial_division_by_higher_degree_returns_zero_quotient() {
    let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5]));
    let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));

    let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should succeed");

    assert!(quotient.is_zero());
    assert_eq!(remainder, dividend);
}

#[test]
fn dense_polynomial_division_handles_exact_division() {
    let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
    let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]));

    let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should succeed");

    assert_eq!(
        quotient,
        DensePolynomial::<F17>::new(f17_coefficients(&[2, 1]))
    );
    assert!(remainder.is_zero());
}

#[test]
fn dense_polynomial_division_returns_expected_remainder() {
    let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 0, 1]));
    let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));

    let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should succeed");

    assert_eq!(
        quotient,
        DensePolynomial::<F17>::new(f17_coefficients(&[0, 1]))
    );
    assert_eq!(
        remainder,
        DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]))
    );
}

#[test]
fn dense_polynomial_division_by_constant_scales_coefficients() {
    let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[4, 8, 12]));
    let divisor = DensePolynomial::<F17>::constant(F17::from_i64(4));

    let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should succeed");

    assert!(remainder.is_zero());
    assert_eq!(
        quotient,
        DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]))
    );
}

#[test]
fn dense_polynomial_division_by_zero_is_rejected() {
    let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]));
    let divisor = DensePolynomial::<F17>::new(Vec::new());

    assert_eq!(
        dividend.div_rem(&divisor),
        Err(PolynomialError::DivisionByZeroPolynomial)
    );
}

#[test]
fn dense_polynomial_gcd_returns_common_factor_in_monic_form() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 3, 3, 1]));

    let gcd = lhs.gcd(&rhs);

    assert_eq!(gcd, DensePolynomial::<F17>::new(f17_coefficients(&[1, 1])));
}

#[test]
fn dense_polynomial_gcd_of_coprime_inputs_is_one() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]));

    let gcd = lhs.gcd(&rhs);

    assert_eq!(gcd, DensePolynomial::<F17>::constant(F17::one()));
}

#[test]
fn dense_polynomial_gcd_handles_zero_inputs() {
    let zero = DensePolynomial::<F17>::new(Vec::new());
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[2, 4]));

    assert_eq!(zero.gcd(&zero), DensePolynomial::<F17>::new(Vec::new()));
    assert_eq!(
        zero.gcd(&polynomial),
        DensePolynomial::<F17>::new(f17_coefficients(&[9, 1]))
    );
    assert_eq!(
        polynomial.gcd(&zero),
        DensePolynomial::<F17>::new(f17_coefficients(&[9, 1]))
    );
}

#[test]
fn dense_polynomial_partial_eq_respects_rational_normalization() {
    let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (1, 3)]));
    let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(2, 4), (3, 9)]));

    assert_eq!(lhs, rhs);
}

#[test]
fn dense_polynomial_mul_by_zero_returns_zero() {
    let lhs = DensePolynomial::<F17>::new(Vec::new());
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]));

    assert!(lhs.mul(&rhs).is_zero());
    assert!(rhs.mul(&lhs).is_zero());
}

#[test]
fn dense_polynomial_add_and_constant_helpers_compose_naturally() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5]));
    let rhs = DensePolynomial::<F17>::constant(F17::from_i64(1));
    let sum = lhs.add(&rhs);

    assert_eq!(sum, DensePolynomial::<F17>::new(f17_coefficients(&[4, 5])));
    assert!(DensePolynomial::<F17>::constant(F17::one()).is_monic());
}

#[test]
fn dense_polynomial_formal_derivative_matches_univariate_trait_call() {
    let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[6, 5, 4]));

    assert_eq!(
        <DensePolynomial<F17> as UnivariatePolynomial<F17>>::derivative(&polynomial),
        DensePolynomial::<F17>::new(f17_coefficients(&[5, 8]))
    );
}

#[test]
fn dense_polynomial_trait_gcd_agrees_with_inherent_gcd() {
    let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
    let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 3, 3, 1]));

    assert_eq!(
        <DensePolynomial<F17> as UnivariatePolynomial<F17>>::gcd(&lhs, &rhs),
        DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]))
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn dense_polynomial_gcd_is_commutative_over_small_f17_inputs(
        lhs in arb_dense_polynomial::<F17>(PolynomialStrategyConfig::default()),
        rhs in arb_dense_polynomial::<F17>(PolynomialStrategyConfig::default()),
    ) {
        prop_assert_eq!(lhs.gcd(&rhs), rhs.gcd(&lhs));
    }

    #[test]
    fn dense_polynomial_mul_agrees_with_sparse_conversion(
        lhs_coefficients in prop::collection::vec(0u64..17, 0..6),
        rhs_coefficients in prop::collection::vec(0u64..17, 0..6),
    ) {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&lhs_coefficients));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&rhs_coefficients));

        let lhs_sparse = SparsePolynomial::<F17>::new(
            lhs.coefficients()
                .iter()
                .cloned()
                .enumerate()
                .map(|(degree, coefficient)| SparsePolynomialTerm { coefficient, degree })
                .collect(),
        );
        let rhs_sparse = SparsePolynomial::<F17>::new(
            rhs.coefficients()
                .iter()
                .cloned()
                .enumerate()
                .map(|(degree, coefficient)| SparsePolynomialTerm { coefficient, degree })
                .collect(),
        );

        let lhs_dense_roundtrip = DensePolynomial::<F17>::from(lhs_sparse);
        let rhs_dense_roundtrip = DensePolynomial::<F17>::from(rhs_sparse);

        prop_assert!(lhs == lhs_dense_roundtrip);
        prop_assert!(rhs == rhs_dense_roundtrip);
        prop_assert!(lhs.mul(&rhs) == lhs_dense_roundtrip.mul(&rhs_dense_roundtrip));
    }

    #[test]
    fn property_dense_and_sparse_evaluation_agree(
        polynomial in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
            max_len: 6,
            ..PolynomialStrategyConfig::default()
        }),
        point in prop::num::u64::ANY.prop_map(F17::from_i64),
    ) {
        let sparse = SparsePolynomial::<F17>::new(
            polynomial
                .coefficients()
                .iter()
                .enumerate()
                .map(|(degree, coefficient)| SparsePolynomialTerm {
                    coefficient: *coefficient,
                    degree,
                })
                .collect(),
        );
        let dense_value = polynomial.evaluate(&point).expect("dense evaluation should succeed");
        let sparse_value = sparse.evaluate(&point).expect("sparse evaluation should succeed");

        prop_assert!(F17::eq(&dense_value, &sparse_value));
    }

    #[test]
    fn property_lagrange_interpolation_recovers_small_dense_polynomials(
        case in interpolation_case(),
    ) {
        let (polynomial, samples) = case;
        let interpolated = DensePolynomial::<F17>::lagrange_interpolate(&samples)
            .expect("interpolation should succeed");
        prop_assert_eq!(interpolated, polynomial);
    }
}
