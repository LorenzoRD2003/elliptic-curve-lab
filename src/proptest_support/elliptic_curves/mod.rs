//! Elliptic-curve-oriented strategies and fixtures.

pub mod analytic;
pub mod division_polynomials;
pub mod endomorphisms;
pub mod frobenius;
pub mod function_fields;
pub mod general_weierstrass;
pub mod points;
pub mod projective;
pub mod short_weierstrass;

pub use analytic::{
    arb_complex_lattice, arb_fundamental_coordinate, arb_interior_fundamental_coordinate,
    arb_stable_real_split_analytic_curve, arb_upper_half_plane_point,
};
pub use division_polynomials::{DivisionPolynomialCase, arb_division_polynomial_case};
pub use endomorphisms::{EndomorphismReportCase, arb_endomorphism_report_case};
pub use frobenius::{FrobeniusCurveCase, arb_frobenius_curve_case};
pub use function_fields::{
    FunctionFieldCase, FunctionFieldPairCase, arb_short_weierstrass_function_case,
    arb_short_weierstrass_function_pair_case,
};
pub use general_weierstrass::{
    arb_general_weierstrass_curve_and_point, arb_nonsingular_general_weierstrass_curve,
};
pub use points::arb_curve_and_point;
pub use projective::{
    arb_general_weierstrass_projective_equivalence_class, arb_general_weierstrass_projective_pair,
    arb_general_weierstrass_projective_point, arb_short_weierstrass_projective_equivalence_class,
    arb_short_weierstrass_projective_pair, arb_short_weierstrass_projective_point,
    rescale_projective_point,
};
pub use short_weierstrass::arb_nonsingular_curve;

pub(crate) fn touch_elliptic_curve_inventory() {
    let analytic_config = crate::proptest_support::config::AnalyticStrategyConfig::default();
    let curve_config = crate::proptest_support::config::CurveStrategyConfig::default();

    let _ = arb_upper_half_plane_point(analytic_config);
    let _ = arb_complex_lattice(analytic_config);
    let _ = arb_fundamental_coordinate();
    let _ = arb_interior_fundamental_coordinate();
    let _ = arb_stable_real_split_analytic_curve();
    let _ = arb_nonsingular_curve::<17>(curve_config);
    let _ = arb_nonsingular_general_weierstrass_curve::<17>(curve_config);
    let _ = arb_curve_and_point::<17>(curve_config);
    let _ = arb_general_weierstrass_curve_and_point::<17>(curve_config);
    let _ = arb_general_weierstrass_projective_point::<17>(curve_config);
    let _ = arb_general_weierstrass_projective_equivalence_class::<17>(curve_config);
    let _ = arb_general_weierstrass_projective_pair::<17>(curve_config);
    let _ = arb_short_weierstrass_projective_point::<17>(curve_config);
    let _ = arb_short_weierstrass_projective_equivalence_class::<17>(curve_config);
    let _ = arb_short_weierstrass_projective_pair::<17>(curve_config);
    let _ = arb_frobenius_curve_case::<17>(curve_config);
    let _ = arb_endomorphism_report_case::<17>(curve_config);
    let _ = arb_division_polynomial_case::<17>(curve_config);
    let _ = arb_short_weierstrass_function_case::<17>(
        curve_config,
        crate::proptest_support::config::PolynomialStrategyConfig::default(),
    );
    let _ = arb_short_weierstrass_function_pair_case::<17>(
        curve_config,
        crate::proptest_support::config::PolynomialStrategyConfig::default(),
    );
    let _ = core::mem::size_of::<FrobeniusCurveCase<17>>();
    let _ = core::mem::size_of::<EndomorphismReportCase<17>>();
    let _ = core::mem::size_of::<DivisionPolynomialCase<17>>();
    let _ = core::mem::size_of::<FunctionFieldCase<17>>();
    let _ = core::mem::size_of::<FunctionFieldPairCase<17>>();
    frobenius::touch_frobenius_case_fields();
    endomorphisms::touch_endomorphism_case_fields();
    division_polynomials::touch_division_polynomial_case_fields();
    function_fields::touch_function_field_case_fields();
}
