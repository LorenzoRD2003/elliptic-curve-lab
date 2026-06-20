use crate::elliptic_curves::short_weierstrass::schoof::ReducedCurveQuotient;
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

/// One reduced class of the form `A(x) + y B(x)` in `F[x,y]/(y^2 - f(x), g(x))`.
///
/// This value type stores one canonical representative of the quotient class
/// `A(x) + y B(x)`, with both `A` and `B` already reduced modulo the chosen
/// univariate modulus `g(x)`.
///
/// Arithmetic is performed relative to one [`ReducedCurveQuotient`]. The
/// quotient provides the extra relation `g(x) = 0`, while the curve provides
/// the short-Weierstrass relation `y^2 = f(x) = x^3 + ax + b`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct ReducedCurveFunction<F: FiniteField> {
    x_part: DensePolynomial<F>,
    y_part: DensePolynomial<F>,
}

impl<F: FiniteField> ReducedCurveFunction<F> {
    #[cfg_attr(not(test), allow(dead_code))]
    fn new_without_reduce(x_part: DensePolynomial<F>, y_part: DensePolynomial<F>) -> Self {
        Self { x_part, y_part }
    }

    /// Builds one canonical reduced representative of `A(x) + y B(x)`.
    ///
    /// Both input polynomials are reduced modulo the active univariate modulus
    /// `g(x)` before storage.
    ///
    /// Complexity: if `m = deg g`, `a = deg A`, and `b = deg B`, then under
    /// the current dense backend this constructor costs `Θ(m(a + b))` field
    /// operations, coming from two remainder reductions modulo `g(x)`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn new(
        quotient: &ReducedCurveQuotient<F>,
        x_part: DensePolynomial<F>,
        y_part: DensePolynomial<F>,
    ) -> Self {
        Self::new_without_reduce(quotient.reduce_poly(&x_part), quotient.reduce_poly(&y_part))
    }

    /// Returns the reduced polynomial `A(x)` in the stored class
    /// `A(x) + y B(x)`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn x_part(&self) -> &DensePolynomial<F> {
        &self.x_part
    }

    /// Returns the reduced polynomial `B(x)` in the stored class
    /// `A(x) + y B(x)`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn y_part(&self) -> &DensePolynomial<F> {
        &self.y_part
    }

    /// Returns whether both reduced components vanish.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn is_zero(&self) -> bool {
        self.x_part.is_zero() && self.y_part.is_zero()
    }

    /// Adds two reduced classes in `F[x, y] / (y^2 - f(x), g(x))`.
    ///
    /// The sum is computed componentwise:
    /// `(A1 + y B1) + (A2 + y B2) = (A1 + A2) + y(B1 + B2)`.
    ///
    /// Because both inputs are already reduced modulo `g(x)`, their component
    /// sums are still valid reduced representatives.
    ///
    /// Complexity: `Θ(m)` field operations on degree-`< m` representatives.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn add(&self, _quotient: &ReducedCurveQuotient<F>, rhs: &Self) -> Self {
        Self::new_without_reduce(self.x_part.add(&rhs.x_part), self.y_part.add(&rhs.y_part))
    }

    /// Negates one reduced class. This sends `A + yB` to `-A + y(-B)`.
    ///
    /// As with addition, negation preserves the degree bounds of the stored
    /// representatives, so no extra reduction modulo `g(x)` is needed.
    ///
    /// Complexity: `Θ(m)` field operations on degree-`< m` representatives.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn neg(&self, _quotient: &ReducedCurveQuotient<F>) -> Self {
        Self::new_without_reduce(self.x_part.neg(), self.y_part.neg())
    }

    /// Subtracts two reduced classes in `F[x, y] / (y^2 - f(x), g(x))`.
    ///
    /// Complexity: `Θ(m)` field operations on degree-`< m` representatives.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn sub(&self, quotient: &ReducedCurveQuotient<F>, rhs: &Self) -> Self {
        self.add(quotient, &rhs.neg(quotient))
    }

    /// Multiplies two reduced classes in `F[x, y] / (y^2 - f(x), g(x))`.
    ///
    /// Writing the factors as `A1(x) + y B1(x)` and `A2(x) + y B2(x)`, their
    /// product is
    ///
    /// `(A1A2 + f B1B2) + y(A1B2 + A2B1)`,
    ///
    /// where `f(x) = x^3 + ax + b` is the short-Weierstrass cubic of the
    /// ambient curve. Both resulting components are then reduced modulo `g(x)`.
    ///
    /// Complexity: a constant number of dense polynomial multiplications and
    /// additions on degree-`< m` representatives, followed by two remainder
    /// reductions modulo `g(x)`. Under the current backend this is `Θ(m^2)`
    /// field operations.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn mul(&self, quotient: &ReducedCurveQuotient<F>, rhs: &Self) -> Self {
        let cubic = quotient.curve().to_cubic();
        let x_part = self
            .x_part
            .mul(&rhs.x_part)
            .add(&cubic.mul(&self.y_part.mul(&rhs.y_part)));
        let y_part = self
            .x_part
            .mul(&rhs.y_part)
            .add(&rhs.x_part.mul(&self.y_part));

        Self::new(quotient, x_part, y_part)
    }
}
