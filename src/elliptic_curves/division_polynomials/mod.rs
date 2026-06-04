//! Educational scaffolding for short-Weierstrass division polynomials and
//! rational torsion helpers.

mod error;
mod polynomial;
mod torsion;

pub use error::DivisionPolynomialError;
pub use polynomial::{
    DivisionPolynomial, DivisionPolynomialForm, DivisionPolynomialXCriterionKind,
    division_polynomial, division_polynomial_base, division_polynomial_x_criterion_kind,
    evaluate_division_polynomial_at_point, evaluate_division_polynomial_x_criterion,
    evaluate_even_division_polynomial_factor_at_x, evaluate_odd_division_polynomial_at_x,
    even_division_polynomial, odd_division_polynomial,
};
pub use torsion::{
    TorsionComparisonReport, compare_division_polynomial_torsion_with_enumeration,
    exact_n_torsion_points_from_division_polynomial, rational_x_candidates_for_division_polynomial,
    torsion_candidates_from_division_polynomial, torsion_points_from_division_polynomial,
};
