use core::fmt;

use crate::fields::traits::Field;

/// General Weierstrass curve model
///
/// `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`.
///
/// This staged initial surface is intentionally small: it records the five
/// coefficients and exposes a human-readable equation string. Validation,
/// invariants, group law, and reduction to short Weierstrass form are added in
/// later milestones.
pub struct GeneralWeierstrassCurve<F: Field> {
    a1: F::Elem,
    a2: F::Elem,
    a3: F::Elem,
    a4: F::Elem,
    a6: F::Elem,
}

impl<F: Field> GeneralWeierstrassCurve<F> {
    /// Builds a general Weierstrass curve descriptor without yet enforcing the
    /// later-stage non-singularity checks.
    pub fn new(a1: F::Elem, a2: F::Elem, a3: F::Elem, a4: F::Elem, a6: F::Elem) -> Self {
        Self { a1, a2, a3, a4, a6 }
    }

    /// Returns the coefficient `a1`.
    pub fn a1(&self) -> &F::Elem {
        &self.a1
    }

    /// Returns the coefficient `a2`.
    pub fn a2(&self) -> &F::Elem {
        &self.a2
    }

    /// Returns the coefficient `a3`.
    pub fn a3(&self) -> &F::Elem {
        &self.a3
    }

    /// Returns the coefficient `a4`.
    pub fn a4(&self) -> &F::Elem {
        &self.a4
    }

    /// Returns the coefficient `a6`.
    pub fn a6(&self) -> &F::Elem {
        &self.a6
    }

    /// Returns the defining equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!(
            "y^2 + ({})xy + ({})y = x^3 + ({})x^2 + ({})x + ({})",
            self.a1, self.a3, self.a2, self.a4, self.a6
        )
    }
}
