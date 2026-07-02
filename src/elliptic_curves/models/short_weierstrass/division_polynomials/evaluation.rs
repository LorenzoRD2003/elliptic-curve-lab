use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::DivisionPolynomialError, traits::CurveModel,
};
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Evaluates `ψ_n(P)` at a finite affine point `P`.
    ///
    /// - if `n` is odd, the implementation evaluates `ψ_n(x)` at the point's
    ///   `x`-coordinate
    /// - if `n` is even, it evaluates the factor `ε_n(x)` and then multiplies by
    ///   the point's `y`-coordinate
    ///
    /// The current implementation does not support evaluating at the point at infinity.
    pub fn evaluate_division_polynomial_at_point(
        &self,
        n: usize,
        point: &AffinePoint<F>,
    ) -> Result<F::Elem, DivisionPolynomialError> {
        if !self.contains(point) {
            return Err(DivisionPolynomialError::Curve(CurveError::PointNotOnCurve));
        }

        let AffinePoint::Finite { x, y } = point else {
            return Err(DivisionPolynomialError::PointAtInfinityNotSupported);
        };

        let criterion = self.evaluate_division_polynomial_x_criterion(n, x)?;

        if self.division_polynomial_uses_y_factor(n)? {
            Ok(F::mul(y, &criterion))
        } else {
            Ok(criterion)
        }
    }

    /// Evaluates an odd division polynomial `ψ_n ∈ F[x]` at an `x`-coordinate.
    /// - `n` must be odd
    /// - `n >= 1`
    pub(crate) fn evaluate_odd_division_polynomial_at_x(
        &self,
        n: usize,
        x: &F::Elem,
    ) -> Result<F::Elem, DivisionPolynomialError> {
        let polynomial = self.odd_division_polynomial(n)?;
        polynomial.evaluate(x).map_err(Into::into)
    }

    /// Evaluates the `F[x]` factor `ε_n(x)` in an even division polynomial
    /// `ψ_n = y ε_n(x)`.
    /// - `n` must be even
    /// - `n >= 2`
    pub(crate) fn evaluate_even_division_polynomial_factor_at_x(
        &self,
        n: usize,
        x: &F::Elem,
    ) -> Result<F::Elem, DivisionPolynomialError> {
        let polynomial = self.even_division_polynomial_factor(n)?;
        polynomial.evaluate(x).map_err(Into::into)
    }

    /// Evaluates the division-polynomial `x`-criterion attached to index `n`.
    ///
    /// This is the right helper when callers only have an `x`-coordinate:
    ///
    /// - if `n` is odd, it evaluates `ψ_n(x)`
    /// - if `n` is even, it evaluates the stripped factor `ε_n(x)` from
    ///   `ψ_n = y ε_n(x)`
    pub(crate) fn evaluate_division_polynomial_x_criterion(
        &self,
        n: usize,
        x: &F::Elem,
    ) -> Result<F::Elem, DivisionPolynomialError> {
        if self.division_polynomial_uses_y_factor(n)? {
            self.evaluate_even_division_polynomial_factor_at_x(n, x)
        } else {
            self.evaluate_odd_division_polynomial_at_x(n, x)
        }
    }
}
