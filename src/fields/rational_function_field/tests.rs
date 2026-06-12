use proptest::prelude::*;

use crate::fields::{
    Field, FieldError, Fp, PthRootExtraction, Q, RationalFunction, RationalFunctionField,
};
use crate::polynomials::DensePolynomial;
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::fields::arb_rational_function;
use crate::proptest_support::polynomials::arb_dense_polynomial;

type F17 = Fp<17>;
type F17x = RationalFunctionField<F17>;
type Qx = RationalFunctionField<Q>;

crate::fields::define_fp_quadratic_extension!(
    spec: F17Sqrt3RationalFunctionFrobeniusSpec,
    field: F17Sqrt3RationalFunctionFrobenius,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3)) for rational-function Frobenius tests",
);

fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
    DensePolynomial::<F17>::new(values.iter().copied().map(F17::elem_from_u64).collect())
}

fn q_dense(values: &[(i64, i64)]) -> DensePolynomial<Q> {
    DensePolynomial::<Q>::new(
        values
            .iter()
            .map(|&(numerator, denominator)| {
                Q::div(&Q::from_i64(numerator), &Q::from_i64(denominator))
                    .expect("denominator should be non-zero")
            })
            .collect(),
    )
}

#[test]
fn rational_function_constructor_rejects_zero_denominator() {
    let error =
        RationalFunction::<F17>::new(f17_dense(&[1, 2]), DensePolynomial::<F17>::new(Vec::new()))
            .expect_err("zero denominator should fail");

    assert_eq!(error, FieldError::DivisionByZero);
}

#[test]
fn rational_function_constructor_reduces_and_makes_denominator_monic() {
    let function = RationalFunction::<F17>::new(f17_dense(&[2, 6, 4]), f17_dense(&[2, 2]))
        .expect("rational function should exist");

    assert_eq!(function.numerator(), &f17_dense(&[1, 2]));
    assert_eq!(function.denominator(), &f17_dense(&[1]));
    assert!(function.denominator().is_monic());
}

#[test]
fn rational_function_zero_is_stored_as_zero_over_one() {
    let function =
        RationalFunction::<F17>::new(DensePolynomial::<F17>::new(Vec::new()), f17_dense(&[3, 4]))
            .expect("zero rational function should exist");

    assert!(function.is_zero());
    assert_eq!(
        function.numerator(),
        &DensePolynomial::<F17>::new(Vec::new())
    );
    assert_eq!(function.denominator(), &f17_dense(&[1]));
}

#[test]
fn rational_function_addition_returns_canonical_result() {
    let left =
        RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1])).expect("left exists");
    let right =
        RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 16])).expect("right exists");

    let sum = left.add(&right);

    assert_eq!(sum.numerator(), &f17_dense(&[15]));
    assert_eq!(sum.denominator(), &f17_dense(&[16, 0, 1]));
}

#[test]
fn rational_function_multiplication_and_inverse_work() {
    let left =
        RationalFunction::<F17>::new(f17_dense(&[1, 1]), f17_dense(&[1, 2])).expect("left exists");
    let right =
        RationalFunction::<F17>::new(f17_dense(&[1, 3]), f17_dense(&[1, 4])).expect("right exists");

    let product = left.mul(&right);
    assert_eq!(product.numerator(), &f17_dense(&[15, 9, 11]));
    assert_eq!(product.denominator(), &f17_dense(&[15, 5, 1]));

    let inverse = left
        .inverse()
        .expect("non-zero rational function is invertible");
    assert_eq!(
        left.mul(&inverse),
        RationalFunction::<F17>::constant(F17::one())
    );
}

#[test]
fn rational_function_division_rejects_zero_divisor() {
    let left =
        RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1])).expect("left exists");
    let zero = RationalFunction::<F17>::constant(F17::zero());

    let error = left
        .div(&zero)
        .expect_err("division by zero rational function should fail");
    assert_eq!(error, FieldError::DivisionByZero);
}

#[test]
fn rational_function_derivative_uses_the_quotient_rule() {
    let function =
        RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1, 1])).expect("exists");

    let derivative = function.derivative();

    assert_eq!(derivative.numerator(), &f17_dense(&[1]));
    assert_eq!(derivative.denominator(), &f17_dense(&[1, 2, 1]));
}

#[test]
fn rational_function_equality_uses_canonical_reduction() {
    let left =
        RationalFunction::<F17>::new(f17_dense(&[2, 2]), f17_dense(&[2])).expect("left exists");
    let right = RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 1]));

    assert_eq!(left, right);
    assert!(right.denominator().is_monic());
}

#[test]
fn rational_function_derivative_works_over_q_too() {
    let function =
        RationalFunction::<Q>::new(q_dense(&[(0, 1), (1, 1)]), q_dense(&[(1, 1), (1, 1)]))
            .expect("exists");

    let derivative = function.derivative();

    assert_eq!(derivative.numerator(), &q_dense(&[(1, 1)]));
    assert_eq!(
        derivative.denominator(),
        &q_dense(&[(1, 1), (2, 1), (1, 1)])
    );
}

#[test]
fn rational_function_pth_root_recovers_a_prime_field_example() {
    let mut numerator = vec![F17::zero(); 18];
    numerator[17] = F17::one();
    let mut denominator = vec![F17::zero(); 18];
    denominator[0] = F17::one();
    denominator[17] = F17::one();

    let function = RationalFunction::<F17>::new(
        DensePolynomial::<F17>::new(numerator),
        DensePolynomial::<F17>::new(denominator),
    )
    .expect("example rational function should exist");

    let root = function
        .pth_root()
        .expect("x^17 / (1 + x^17) should be a 17-th power in F17(x)");

    assert_eq!(
        root,
        RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1, 1]))
            .expect("x / (1 + x) should exist")
    );
    assert!(function.has_pth_root());
}

#[test]
fn rational_function_pth_root_rejects_non_pth_power_numerator() {
    let function = RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1]))
        .expect("x should define a rational function");

    assert_eq!(function.pth_root(), None);
    assert!(!function.has_pth_root());
}

#[test]
fn rational_function_pth_root_rejects_non_pth_power_denominator() {
    let function = RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1]))
        .expect("1 / (1 + x) should define a rational function");

    assert_eq!(function.pth_root(), None);
    assert!(!function.has_pth_root());
}

#[test]
fn rational_function_zero_has_a_pth_root() {
    let zero = RationalFunction::<F17>::constant(F17::zero());

    assert_eq!(zero.pth_root(), Some(zero.clone()));
    assert!(zero.has_pth_root());
}

#[test]
fn rational_function_inverts_absolute_frobenius_pullback_in_prime_field_case() {
    let function = RationalFunction::<F17>::new(
        {
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[0] = F17::elem_from_u64(3);
            coefficients[17] = F17::elem_from_u64(2);
            DensePolynomial::<F17>::new(coefficients)
        },
        {
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[0] = F17::one();
            coefficients[17] = F17::elem_from_u64(4);
            DensePolynomial::<F17>::new(coefficients)
        },
    )
    .expect("example should exist");

    let preimage = function
        .inverse_absolute_frobenius_pullback_from_twist()
        .expect("x -> x^p image should invert");

    assert_eq!(
        preimage,
        RationalFunction::<F17>::new(f17_dense(&[3, 2]), f17_dense(&[1, 4]))
            .expect("preimage should exist")
    );
}

#[test]
fn rational_function_inverts_absolute_frobenius_pullback_over_extension_field() {
    let alpha = F17Sqrt3RationalFunctionFrobenius::element(vec![F17::zero(), F17::one()]);
    let function = RationalFunction::<F17Sqrt3RationalFunctionFrobenius>::new(
        {
            let mut coefficients = vec![F17Sqrt3RationalFunctionFrobenius::zero(); 18];
            coefficients[0] = alpha.clone();
            coefficients[17] = F17Sqrt3RationalFunctionFrobenius::one();
            DensePolynomial::<F17Sqrt3RationalFunctionFrobenius>::new(coefficients)
        },
        {
            let mut coefficients = vec![F17Sqrt3RationalFunctionFrobenius::zero(); 18];
            coefficients[0] = F17Sqrt3RationalFunctionFrobenius::one();
            coefficients[17] = alpha.clone();
            DensePolynomial::<F17Sqrt3RationalFunctionFrobenius>::new(coefficients)
        },
    )
    .expect("example should exist");

    let preimage = function
        .inverse_absolute_frobenius_pullback_from_twist()
        .expect("extension-field image should invert");

    assert_eq!(
        preimage,
        RationalFunction::<F17Sqrt3RationalFunctionFrobenius>::new(
            DensePolynomial::<F17Sqrt3RationalFunctionFrobenius>::new(vec![
                alpha.clone(),
                F17Sqrt3RationalFunctionFrobenius::one(),
            ]),
            DensePolynomial::<F17Sqrt3RationalFunctionFrobenius>::new(vec![
                F17Sqrt3RationalFunctionFrobenius::one(),
                alpha,
            ]),
        )
        .expect("preimage should exist")
    );
}

#[test]
fn rational_function_inverse_absolute_frobenius_pullback_rejects_non_multiple_degree_terms() {
    let function = RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1]))
        .expect("x should define a rational function");

    assert_eq!(
        function.inverse_absolute_frobenius_pullback_from_twist(),
        None
    );
}

#[test]
fn rational_function_field_zero_one_and_constants_match_value_layer() {
    assert_eq!(F17x::zero(), RationalFunction::<F17>::constant(F17::zero()));
    assert_eq!(F17x::one(), RationalFunction::<F17>::constant(F17::one()));
    assert_eq!(
        F17x::elem_from_u64(5),
        RationalFunction::<F17>::constant(F17::elem_from_u64(5))
    );
}

#[test]
fn rational_function_field_indeterminate_is_x_over_one() {
    let x = F17x::indeterminate();

    assert_eq!(x.numerator(), &f17_dense(&[0, 1]));
    assert_eq!(x.denominator(), &f17_dense(&[1]));
}

#[test]
fn rational_function_field_inverse_rejects_zero() {
    assert_eq!(
        F17x::inverse(&F17x::zero()),
        Err(FieldError::DivisionByZero)
    );
    assert_eq!(F17x::inv(&F17x::zero()), None);
}

#[test]
fn rational_function_field_inherits_characteristic_and_is_not_algebraically_closed() {
    const _: () = {
        assert!(!F17x::IS_ALGEBRAICALLY_CLOSED);
        assert!(!Qx::IS_ALGEBRAICALLY_CLOSED);
    };

    assert_eq!(F17x::characteristic(), F17::characteristic());
    assert_eq!(Qx::characteristic(), Q::characteristic());
}

#[test]
fn rational_function_field_operations_delegate_to_rational_function_values() {
    let x = F17x::indeterminate();
    let one = F17x::one();
    let one_over_x = F17x::inverse(&x).expect("x should be invertible in F(x)");

    assert_eq!(F17x::mul(&one_over_x, &x), one);
    assert_eq!(F17x::add(&x, &one), x.add(&one));
    assert_eq!(F17x::sub(&x, &one), x.sub(&one));
    assert_eq!(F17x::neg(&x), x.neg());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(40))]

    #[test]
    fn rational_function_canonicalization_forgets_common_polynomial_factors(
        function in arb_rational_function::<F17>(PolynomialStrategyConfig {
            max_len: 5,
            ..PolynomialStrategyConfig::default()
        }),
        common_factor in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
            max_len: 4,
            require_nonzero_leading_coefficient: true,
            ..PolynomialStrategyConfig::default()
        }),
    ) {
        prop_assume!(!common_factor.is_zero());

        let redundant_numerator = function.numerator().mul(&common_factor);
        let redundant_denominator = function.denominator().mul(&common_factor);

        let reconstructed = RationalFunction::<F17>::new(redundant_numerator, redundant_denominator)
            .expect("multiplying numerator and denominator by the same non-zero polynomial should keep a valid presentation");

        prop_assert_eq!(reconstructed.clone(), function);
        prop_assert!(reconstructed.denominator().is_monic());
    }
}
