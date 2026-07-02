use crate::fields::traits::*;
use crate::polynomials::DensePolynomial;

/// For a short-Weierstrass curve `E: y^2 = x^3 + ax + b`,
/// the first division polynomials are:
///
/// - `ψ_0 = 0`
/// - `ψ_1 = 1`
/// - `ψ_2 = 2y`
/// - `ψ_3 = 3x^4 + 6ax^2 + 12bx - a^2`
/// - `ψ_4 = 4y (x^6 + 5ax^4 + 20bx^3 - 5a^2x^2 - 4abx - 8b^2 - a^3)`
///
/// The important structural distinction is:
/// - for odd `n`, `ψ_n ∈ F[x]`
/// - for even `n`, `ψ_n ∈ yF[x]`
///
/// This enum therefore records whether the object lives directly in `F[x]`
/// or in `yF[x]`, while reusing `DensePolynomial<F>` for the polynomial part
/// in the `x`-coordinate.
#[derive(Clone, Debug)]
pub enum DivisionPolynomialForm<F: Field> {
    /// A division polynomial that lies in `F[x]`.
    InX(DensePolynomial<F>),
    /// A division polynomial of the form `y * f(x)` with `f(x) ∈ F[x]`.
    YTimes(DensePolynomial<F>),
}

impl<F: Field> PartialEq for DivisionPolynomialForm<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::InX(lhs), Self::InX(rhs)) => lhs == rhs,
            (Self::YTimes(lhs), Self::YTimes(rhs)) => lhs == rhs,
            _ => false,
        }
    }
}

impl<F: Field> DivisionPolynomialForm<F> {
    /// Returns whether this division polynomial lies directly in `F[x]`.
    pub fn is_x_polynomial(&self) -> bool {
        matches!(self, Self::InX(_))
    }

    /// Returns whether this division polynomial has the form `y * f(x)`.
    pub fn is_y_times_x_polynomial(&self) -> bool {
        matches!(self, Self::YTimes(_))
    }

    /// Returns the stored polynomial factor in `F[x]`.
    ///
    /// For `InX(f)`, this returns `f`.
    ///
    /// For `YTimes(f)`, this returns the factor `f` in the
    /// decomposition `y * f(x)`.
    pub fn x_factor(&self) -> &DensePolynomial<F> {
        match self {
            Self::InX(polynomial) | Self::YTimes(polynomial) => polynomial,
        }
    }

    /// Returns the stored polynomial factor in `F[x]`.
    ///
    /// This is a semantic alias for [`Self::x_factor`] that reads more
    /// naturally when callers only care about the underlying stored
    /// `x`-polynomial and not about whether the honest division-polynomial
    /// shape is `f(x)` or `y f(x)`.
    pub fn stored_x_factor(&self) -> &DensePolynomial<F> {
        self.x_factor()
    }

    /// Wraps a polynomial that lies directly in `F[x]`.
    pub(crate) fn x_polynomial(polynomial: DensePolynomial<F>) -> Self {
        Self::InX(polynomial)
    }

    /// Wraps a polynomial factor `f(x)` to represent `y * f(x)`.
    pub(crate) fn y_times_x_polynomial(polynomial: DensePolynomial<F>) -> Self {
        Self::YTimes(polynomial)
    }
}
