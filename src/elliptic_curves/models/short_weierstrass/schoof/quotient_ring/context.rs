use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, short_weierstrass::schoof::ReducedCurveFunction,
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

/// The reduced short-Weierstrass quotient `F[x, y] / (y^2 - f(x), g(x))`.
///
/// For a short-Weierstrass curve `E: y^2 = f(x)` and a non-zero univariate
/// polynomial `g(x)`, this context packages the extra quotient relation
/// `g(x) = 0` that Schoof-style torsion computations need.
///
/// The module currently provides only the first layer of quotient functionality:
///
/// - reduction of univariate polynomials modulo `g(x)`
/// - canonical representatives of `0`, `1`, `x`, and `y`
/// - partial invertibility information in `F[x] / (g(x))`
///
/// It intentionally does **not** yet provide division or endomorphism
/// arithmetic. Those belong to later milestones.
#[derive(Clone, Debug)]
pub(crate) struct ReducedCurveQuotient<F: FiniteField> {
    #[cfg_attr(not(test), allow(dead_code))]
    curve: ShortWeierstrassCurve<F>,
    modulus: DensePolynomial<F>,
}

impl<F: FiniteField> ReducedCurveQuotient<F> {
    /// Builds the reduced quotient `F[x, y] / (y^2 - f(x), g(x))`.
    ///
    /// The supplied `modulus` is normalized to its monic representative so the
    /// quotient keeps one canonical stored modulus. The zero polynomial is
    /// rejected, because quotienting by `g(x) = 0` would not define the
    /// intended extra relation.
    ///
    /// Complexity: `Θ(deg g)` plus one monic normalization.
    pub(crate) fn new(
        curve: ShortWeierstrassCurve<F>,
        modulus: DensePolynomial<F>,
    ) -> Result<Self, CurveError> {
        if modulus.is_zero() {
            return Err(CurveError::ZeroReducedCurveQuotientModulus);
        }
        let modulus = modulus
            .make_monic()
            .expect("a non-zero polynomial over a field admits monic normalization");

        Ok(Self { curve, modulus })
    }

    /// Returns the ambient short-Weierstrass curve.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.curve
    }

    /// Returns the canonical monic representative of the stored modulus
    /// `g(x)`.
    pub(crate) fn modulus(&self) -> &DensePolynomial<F> {
        &self.modulus
    }

    /// Reduces one univariate polynomial modulo the stored relation `g(x)`.
    ///
    /// The result is the unique remainder of degree strictly less than
    /// `deg g`, or the zero polynomial when the input is divisible by `g`.
    ///
    /// Complexity: one dense Euclidean remainder computation. Under the
    /// current backend this is `Θ(nm)` field operations for input degree `n`
    /// and modulus degree `m`.
    pub(crate) fn reduce_poly(&self, poly: &DensePolynomial<F>) -> DensePolynomial<F> {
        poly.rem(&self.modulus)
            .expect("the stored reduced-curve modulus is non-zero")
    }

    /// Returns the zero class.
    ///
    /// Complexity: `Θ(m)` coefficient work from reducing two zero
    /// polynomials modulo a degree-`m` modulus.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn zero(&self) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            self,
            DensePolynomial::new(Vec::new()),
            DensePolynomial::new(Vec::new()),
        )
    }

    /// Returns the constant-one class.
    ///
    /// Complexity: `Θ(m)` coefficient work.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn one(&self) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            self,
            DensePolynomial::constant(F::one()),
            DensePolynomial::new(Vec::new()),
        )
    }

    /// Returns the distinguished `x` class.
    ///
    /// Complexity: `Θ(m)` coefficient work.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn x(&self) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            self,
            DensePolynomial::new(vec![F::zero(), F::one()]),
            DensePolynomial::new(Vec::new()),
        )
    }

    /// Returns the distinguished `y` class.
    ///
    /// Complexity: `Θ(m)`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn y(&self) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            self,
            DensePolynomial::new(Vec::new()),
            DensePolynomial::constant(F::one()),
        )
    }
}
