use super::{GeneralWeierstrassYFiberError, solve_in_characteristic_two};
use crate::elliptic_curves::{
    GeneralWeierstrassCurve,
    traits::{AffineCurveModel, CurveModel},
};
use crate::fields::{
    Fp,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
    traits::Field,
};

type F2 = Fp<2>;
type F4 = ExtensionField<F4GeneralWeierstrassSpec>;
type F5 = Fp<5>;

#[derive(Clone, Copy)]
struct F4GeneralWeierstrassSpec;

impl ExtensionFieldSpec for F4GeneralWeierstrassSpec {
    type Base = F2;

    const NAME: &'static str = "F4 for general Weierstrass y-fiber tests";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![F2::one(), F2::one(), F2::one()])
            .expect("x^2 + x + 1 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

#[test]
fn y_fiber_equation_records_the_expected_coefficients() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let equation = curve.y_fiber_equation(&F5::zero());

    assert!(F5::eq(equation.u(), &F5::one()));
    assert!(F5::eq(equation.v(), &F5::zero()));
}

#[test]
fn y_fiber_odd_characteristic_completed_square_data_matches_the_expected_values() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let equation = curve.y_fiber_equation(&F5::zero());

    assert!(F5::eq(
        &equation
            .odd_characteristic_shift()
            .expect("characteristic five should support completing the square"),
        &F5::from_i64(3)
    ));
    assert!(F5::eq(
        &equation
            .odd_characteristic_completed_rhs()
            .expect("characteristic five should support completing the square"),
        &F5::from_i64(4)
    ));
}

#[test]
fn y_fiber_solve_in_odd_characteristic_finds_the_two_roots_when_they_exist() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let x = F5::zero();
    let equation = curve.y_fiber_equation(&x);

    let (left_y, right_y) = equation
        .solve_in_odd_characteristic()
        .expect("characteristic five should support completing the square")
        .expect("x = 0 should admit two y-solutions");

    assert!(
        curve.contains(
            &curve
                .point(x, left_y)
                .expect("left root should define a point")
        )
    );
    assert!(
        curve.contains(
            &curve
                .point(x, right_y)
                .expect("right root should define a point")
        )
    );
    assert!(
        (F5::eq(&left_y, &F5::zero()) && F5::eq(&right_y, &F5::from_i64(4)))
            || (F5::eq(&left_y, &F5::from_i64(4)) && F5::eq(&right_y, &F5::zero()))
    );
}

#[test]
fn y_fiber_solve_in_odd_characteristic_detects_a_repeated_root() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let equation = curve.y_fiber_equation(&F5::from_i64(2));

    let (left_y, right_y) = equation
        .solve_in_odd_characteristic()
        .expect("characteristic five should support completing the square")
        .expect("x = 2 should admit one repeated root");

    assert_eq!(left_y, right_y);
    assert!(F5::eq(&left_y, &F5::one()));
}

#[test]
fn y_fiber_solve_in_odd_characteristic_returns_none_when_no_square_root_exists() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let equation = curve.y_fiber_equation(&F5::from_i64(3));

    assert_eq!(
        equation
            .solve_in_odd_characteristic()
            .expect("characteristic five should support completing the square"),
        None
    );
}

#[test]
fn y_fiber_odd_characteristic_formula_rejects_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let equation = curve.y_fiber_equation(&F2::one());

    assert_eq!(
        equation.odd_characteristic_shift(),
        Err(GeneralWeierstrassYFiberError::UnsupportedCharacteristic { characteristic: 2 })
    );
    assert_eq!(
        equation.solve_in_odd_characteristic(),
        Err(GeneralWeierstrassYFiberError::UnsupportedCharacteristic { characteristic: 2 })
    );
}

#[test]
fn y_fiber_characteristic_two_pure_square_case_returns_the_unique_repeated_root() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let equation = curve.y_fiber_equation(&F2::one());

    assert!(F2::is_zero(equation.u()));
    assert!(F2::is_zero(equation.v()));
    assert_eq!(
        solve_in_characteristic_two(&equation)
            .expect("characteristic two should support the pure-square branch"),
        Some((F2::zero(), F2::zero()))
    );
}

#[test]
fn y_fiber_characteristic_two_normalized_rhs_matches_the_expected_value() {
    let curve = GeneralWeierstrassCurve::<F4>::new(
        F4::zero(),
        F4::zero(),
        F4::one(),
        F4::zero(),
        F4::zero(),
    )
    .expect("non-singular curve in characteristic two");
    let equation = curve.y_fiber_equation(&F4::one());

    assert!(!F4::is_zero(equation.u()));
    assert!(F4::eq(
        &equation
            .characteristic_two_normalized_rhs()
            .expect("characteristic two should support Artin-Schreier normalization")
            .expect("u != 0 should produce a normalized rhs"),
        &F4::one()
    ));
}

#[test]
fn y_fiber_solve_in_characteristic_two_returns_none_when_the_artin_schreier_equation_has_no_solution()
 {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let equation = curve.y_fiber_equation(&F2::zero());

    assert_eq!(
        solve_in_characteristic_two(&equation)
            .expect("characteristic two should support the Artin-Schreier branch"),
        None
    );
}

#[test]
fn y_fiber_solve_in_characteristic_two_returns_two_roots_when_the_artin_schreier_equation_is_solvable()
 {
    let curve = GeneralWeierstrassCurve::<F4>::new(
        F4::zero(),
        F4::zero(),
        F4::one(),
        F4::zero(),
        F4::zero(),
    )
    .expect("non-singular curve in characteristic two");
    let equation = curve.y_fiber_equation(&F4::one());

    let (left_y, right_y) = equation
        .solve()
        .expect("characteristic two should support the Artin-Schreier branch")
        .expect("the normalized Artin-Schreier equation should be solvable");

    assert_ne!(left_y, right_y);
    assert!(
        curve.contains(
            &curve
                .point(F4::one(), left_y.clone())
                .expect("left root should define a point")
        )
    );
    assert!(
        curve.contains(
            &curve
                .point(F4::one(), right_y.clone())
                .expect("right root should define a point")
        )
    );
}

#[test]
fn y_fiber_unified_solver_matches_the_odd_characteristic_formula() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let equation = curve.y_fiber_equation(&F5::zero());

    assert_eq!(
        equation.solve().expect("unified solver should succeed"),
        equation
            .solve_in_odd_characteristic()
            .expect("odd-characteristic formula should succeed")
    );
}

#[test]
fn y_fiber_curve_helper_matches_the_equation_solver() {
    let curve = GeneralWeierstrassCurve::<F4>::new(
        F4::zero(),
        F4::zero(),
        F4::one(),
        F4::zero(),
        F4::zero(),
    )
    .expect("non-singular curve in characteristic two");
    let x = F4::one();

    assert_eq!(
        curve
            .solve_y_fiber(&x)
            .expect("curve helper should solve the fiber"),
        curve
            .y_fiber_equation(&x)
            .solve()
            .expect("equation solver should agree")
    );
}
