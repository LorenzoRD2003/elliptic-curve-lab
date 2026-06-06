use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::One;

use crate::fields::{ComplexApprox, Field, FieldError, Fp, Q};
use crate::polynomials::{
    DensePolynomial, IrreducibilityStatus, PolynomialError, ReducibilityReason,
};

use crate::polynomials::irreducibility::{irreducibility_status, is_irreducible};

type F17 = Fp<17>;
type F15 = Fp<15>;

fn dense(values: &[u64]) -> DensePolynomial<F17> {
    DensePolynomial::new(values.iter().copied().map(F17::elem_from_u64).collect())
}

fn complex_dense(values: &[(f64, f64)]) -> DensePolynomial<ComplexApprox> {
    DensePolynomial::new(
        values
            .iter()
            .copied()
            .map(|(re, im)| Complex64::new(re, im))
            .collect(),
    )
}

#[test]
fn constants_are_classified_as_constants() {
    assert_eq!(
        irreducibility_status(&DensePolynomial::<F17>::new(Vec::new()))
            .expect("zero polynomial should classify"),
        IrreducibilityStatus::Constant
    );
    assert_eq!(
        irreducibility_status(&DensePolynomial::<F17>::constant(F17::elem_from_u64(9)))
            .expect("constant polynomial should classify"),
        IrreducibilityStatus::Constant
    );
    assert!(
        !is_irreducible(&DensePolynomial::<F17>::constant(F17::elem_from_u64(9)))
            .expect("constant polynomial should classify")
    );
}

#[test]
fn linear_polynomials_are_classified_as_linear_and_irreducible() {
    let polynomial = dense(&[3, 5]);

    assert_eq!(
        irreducibility_status(&polynomial).expect("linear polynomial should classify"),
        IrreducibilityStatus::Linear
    );
    assert!(is_irreducible(&polynomial).expect("linear polynomial should classify"));
}

#[test]
fn reducible_quadratic_returns_a_witness_factorization() {
    let polynomial = dense(&[1, 0, 1]);
    let status = irreducibility_status(&polynomial).expect("quadratic polynomial should classify");

    match status {
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            assert_eq!(polynomial, divisor.mul(&quotient));
            assert_eq!(divisor.degree(), Some(1));
            assert_eq!(quotient.degree(), Some(1));
        }
        other => panic!("expected reducible status, got {other:?}"),
    }

    assert!(!is_irreducible(&polynomial).expect("quadratic polynomial should classify"));
}

#[test]
fn irreducible_quadratic_is_reported_as_irreducible() {
    let polynomial = dense(&[1, 0, 3]);

    assert_eq!(
        irreducibility_status(&polynomial).expect("quadratic polynomial should classify"),
        IrreducibilityStatus::Irreducible
    );
    assert!(is_irreducible(&polynomial).expect("quadratic polynomial should classify"));
}

#[test]
fn reducibility_witness_is_reported_in_the_original_scaling() {
    let polynomial = dense(&[2, 4, 2]);
    let status = irreducibility_status(&polynomial).expect("quadratic polynomial should classify");

    match status {
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            assert_eq!(polynomial, divisor.mul(&quotient));
        }
        other => panic!("expected reducible status, got {other:?}"),
    }
}

#[test]
fn invalid_prime_field_modulus_is_reported_as_a_typed_error() {
    let polynomial = DensePolynomial::<F15>::new(vec![F15::one(), F15::zero(), F15::one()]);

    assert_eq!(
        irreducibility_status(&polynomial),
        Err(PolynomialError::InvalidBaseField(
            FieldError::InvalidModulus { modulus: 15 }
        ))
    );
}

#[test]
fn complex_constants_and_linears_keep_the_basic_classification_conventions() {
    assert_eq!(
        irreducibility_status(&DensePolynomial::<ComplexApprox>::new(Vec::new()))
            .expect("zero polynomial should classify"),
        IrreducibilityStatus::Constant
    );
    assert_eq!(
        irreducibility_status(&complex_dense(&[(2.0, 0.0), (1.0, -1.0)]))
            .expect("linear complex polynomial should classify"),
        IrreducibilityStatus::Linear
    );
}

#[test]
fn complex_higher_degree_polynomials_are_reducible_by_field_property() {
    let polynomial = complex_dense(&[(1.0, 0.0), (0.0, 0.0), (1.0, 0.0)]);

    assert_eq!(
        irreducibility_status(&polynomial).expect("complex polynomial should classify"),
        IrreducibilityStatus::ReducibleWithoutWitness {
            reason: ReducibilityReason::AlgebraicallyClosedField,
        }
    );
    assert!(!is_irreducible(&polynomial).expect("complex polynomial should classify"));
}

#[test]
fn rationals_can_be_certified_irreducible_on_small_examples() {
    let polynomial =
        DensePolynomial::<Q>::new(vec![Q::from_i64(1), Q::from_i64(0), Q::from_i64(1)]);

    assert_eq!(
        irreducibility_status(&polynomial).expect("x^2 + 1 should be irreducible over Q"),
        IrreducibilityStatus::Irreducible
    );
}

#[test]
fn rational_reducibility_can_be_detected_via_linear_root_search() {
    let polynomial =
        DensePolynomial::<Q>::new(vec![Q::from_i64(-1), Q::from_i64(0), Q::from_i64(1)]);

    match irreducibility_status(&polynomial).expect("x^2 - 1 should be reducible over Q") {
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            assert_eq!(polynomial, divisor.mul(&quotient));
            assert_eq!(divisor.degree(), Some(1));
        }
        other => panic!("expected reducible status, got {other:?}"),
    }
}

#[test]
fn rational_irreducibility_handles_small_denominator_examples_exactly() {
    let half = BigRational::new(BigInt::from(1), BigInt::from(2));
    let polynomial = DensePolynomial::<Q>::new(vec![half.clone(), Q::one(), half]);

    match irreducibility_status(&polynomial).expect("scaled square should be reducible over Q") {
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            assert_eq!(polynomial, divisor.mul(&quotient));
        }
        other => panic!("expected reducible status, got {other:?}"),
    }
}

#[test]
fn rational_backend_reports_inconclusive_cases_honestly() {
    let leading = [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        .into_iter()
        .fold(BigInt::one(), |accumulator, prime| {
            accumulator * BigInt::from(prime)
        });
    let constant = &leading + BigInt::one();
    let polynomial = DensePolynomial::<Q>::new(vec![
        BigRational::from_integer(constant),
        Q::zero(),
        Q::zero(),
        Q::zero(),
        BigRational::from_integer(leading),
    ]);

    assert!(matches!(
        irreducibility_status(&polynomial),
        Err(PolynomialError::UndeterminedIrreducibility(_))
    ));
}
