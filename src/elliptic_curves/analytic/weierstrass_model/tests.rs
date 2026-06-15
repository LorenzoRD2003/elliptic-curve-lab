use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, AnalyticInvariants, AnalyticWeierstrassCurve,
    ApproxTolerance, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint,
};
use crate::fields::{complex_approx::ComplexApprox, traits::Field};

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

#[test]
fn constructor_rejects_nearly_singular_coefficients() {
    assert_eq!(
        AnalyticWeierstrassCurve::new(c(0.0, 0.0), c(0.0, 0.0)),
        Err(AnalyticCurveError::NearlySingularAnalyticCurve)
    );
}

#[test]
fn discriminant_matches_the_analytic_formula() {
    let curve = AnalyticWeierstrassCurve::new(c(12.0, 1.0), c(-3.0, 2.0)).unwrap();
    let expected = AnalyticInvariants::discriminant_from_g2_g3(curve.g2(), curve.g3());

    assert!(ComplexApprox::eq(&curve.discriminant(), &expected));
}

#[test]
fn j_invariant_matches_the_shared_analytic_helper() {
    let curve = AnalyticWeierstrassCurve::new(c(12.0, 1.0), c(4.0, -2.0)).unwrap();
    let expected = AnalyticInvariants::j_invariant_from_g2_g3(curve.g2(), curve.g3()).unwrap();

    assert!(ComplexApprox::eq(&curve.j_invariant().unwrap(), &expected));
}

#[test]
fn from_tau_matches_from_lattice() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let from_tau = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let from_lattice =
        AnalyticWeierstrassCurve::from_lattice(&ComplexLattice::from_tau(tau), truncation).unwrap();

    assert_eq!(from_tau, from_lattice);
}

#[test]
fn membership_report_accepts_curve_points_and_infinity() {
    let curve = AnalyticWeierstrassCurve::new(c(8.0, 0.0), c(-4.0, 0.0)).unwrap();
    let point = AnalyticCurvePoint::new(c(1.0, 0.0), c(0.0, 0.0));
    let infinity = AnalyticCurvePoint::infinity();

    assert!(
        curve
            .membership_report(&point, ApproxTolerance::strict())
            .is_on_curve()
    );
    assert!(
        curve
            .membership_report(&infinity, ApproxTolerance::strict())
            .is_on_curve()
    );
}

#[test]
fn equation_string_mentions_both_invariants() {
    let curve = AnalyticWeierstrassCurve::new(c(8.0, 0.0), c(-4.0, 0.0)).unwrap();
    let equation = curve.equation_string();

    assert!(equation.contains("4x^3"));
    assert!(equation.contains("8"));
    assert!(equation.contains("-4"));
}
