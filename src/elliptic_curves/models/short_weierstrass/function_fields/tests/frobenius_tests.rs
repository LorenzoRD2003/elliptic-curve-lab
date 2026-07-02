use crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunction;
use crate::fields::traits::*;
use crate::fields::{Q, rational_function_field::RationalFunction, traits::PthRootExtraction};

use super::shared::{F17, f17_curve, f17_dense, q_curve};

#[test]
fn derivative_of_x_is_one() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let derivative = field.x().derivative();

    assert!(derivative.is_one());
}

#[test]
fn derivative_of_y_matches_f_prime_over_two_y_in_basis_form() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            Q,
        >::new(q_curve());
    let derivative = field.y().derivative();

    assert!(derivative.a_part().is_zero());
    assert_eq!(
        derivative.b_part(),
        &RationalFunction::<Q>::new(
            super::shared::q_dense(&[(-1, 1), (0, 1), (3, 1)]),
            super::shared::q_dense(&[(0, 1), (-2, 1), (0, 1), (2, 1)])
        )
        .expect("rational function should exist")
    );
}

#[test]
fn derivative_satisfies_product_rule_on_small_example() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let x = field.x();
    let y = field.y();
    let product = x.mul(&y).expect("same-curve multiplication should work");
    let left = product.derivative();
    let right = x
        .derivative()
        .mul(&y)
        .expect("same-curve multiplication should work")
        .add(
            &x.mul(&y.derivative())
                .expect("same-curve multiplication should work"),
        )
        .expect("same-curve addition should work");

    assert_eq!(left, right);
}

#[test]
fn short_weierstrass_function_pth_root_recovers_the_y_generator_from_y_to_the_p() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let y = field.y();
    let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));
    let y_to_the_p = ShortWeierstrassFunction::<F17>::new(
        field.curve().clone(),
        RationalFunction::<F17>::constant(F17::zero()),
        {
            let mut result = RationalFunction::<F17>::constant(F17::one());
            for _ in 0..8 {
                result = result.mul(&rhs);
            }
            result
        },
    );

    assert_eq!(y_to_the_p.pth_root(), Some(y));
    assert!(y_to_the_p.has_pth_root());
}

#[test]
fn short_weierstrass_function_pth_root_rejects_the_x_generator() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());

    assert_eq!(field.x().pth_root(), None);
    assert!(!field.x().has_pth_root());
}

#[test]
fn short_weierstrass_function_pth_root_handles_rational_a_part_examples() {
    let curve = f17_curve();
    let function = ShortWeierstrassFunction::<F17>::from_rational_function(
        curve.clone(),
        RationalFunction::<F17>::new(
            crate::polynomials::DensePolynomial::<F17>::new({
                let mut coefficients = vec![F17::zero(); 18];
                coefficients[17] = F17::one();
                coefficients
            }),
            crate::polynomials::DensePolynomial::<F17>::new({
                let mut coefficients = vec![F17::zero(); 18];
                coefficients[0] = F17::one();
                coefficients[17] = F17::one();
                coefficients
            }),
        )
        .expect("example rational function should exist"),
    );

    let expected = ShortWeierstrassFunction::<F17>::from_rational_function(
        curve,
        RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1, 1]))
            .expect("x / (1 + x) should exist"),
    );

    assert_eq!(function.pth_root(), Some(expected));
    assert!(function.has_pth_root());
}

#[test]
fn short_weierstrass_function_pth_root_recovers_a_mixed_example() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let root = field
        .x()
        .add(&field.y())
        .expect("same-curve addition should work");
    let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));
    let function = ShortWeierstrassFunction::<F17>::new(
        field.curve().clone(),
        RationalFunction::<F17>::from_polynomial({
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[17] = F17::one();
            crate::polynomials::DensePolynomial::<F17>::new(coefficients)
        }),
        {
            let mut result = RationalFunction::<F17>::constant(F17::one());
            for _ in 0..8 {
                result = result.mul(&rhs);
            }
            result
        },
    );

    assert_eq!(function.pth_root(), Some(root));
}
