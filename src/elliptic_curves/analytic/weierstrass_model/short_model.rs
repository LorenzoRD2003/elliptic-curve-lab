use core::fmt;
use num_complex::Complex64;

use crate::elliptic_curves::{
    analytic::AnalyticWeierstrassCurve, short_weierstrass::ShortWeierstrassCurve,
};
use crate::fields::complex_approx::ComplexApprox;

/// Small internal value object for the short-Weierstrass companion
/// `y² = x³ + ax + b` attached to an analytic cubic.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AnalyticShortWeierstrassModel {
    a: Complex64,
    b: Complex64,
    j_invariant: Complex64,
}

impl AnalyticShortWeierstrassModel {
    #[cfg(test)]
    pub(crate) fn from_analytic_curve(curve: &AnalyticWeierstrassCurve) -> Self {
        let short_curve: ShortWeierstrassCurve<ComplexApprox> = curve.as_short_weierstrass();

        Self {
            a: *short_curve.a(),
            b: *short_curve.b(),
            j_invariant: short_curve.j_invariant(),
        }
    }

    pub(crate) fn j_invariant(&self) -> &Complex64 {
        &self.j_invariant
    }
}

impl fmt::Display for AnalyticShortWeierstrassModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "y^2 = x^3 + ({:?})x + ({:?})", self.a, self.b)
    }
}
