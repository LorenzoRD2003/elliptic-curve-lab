use core::fmt;

use crate::fields::{
    FieldError,
    traits::{Field, FiniteField, PthRootExtraction},
};
use crate::polynomials::{DensePolynomial, PolynomialError};

/// Univariate rational function over a field `F`.
///
/// This type models an element of the rational function field `F(x)` through a
/// numerator and denominator in `F[x]`.
///
/// The stored representation is canonicalized eagerly:
///
/// - the denominator must be non-zero
/// - numerator and denominator are divided by their monic gcd
/// - the denominator is normalized to be monic
/// - the zero function is stored as `0 / 1`
///
/// Because of that normalization policy, equality can be understood as
/// equality of rational functions rather than merely equality of one chosen
/// presentation.
pub struct RationalFunction<F: Field> {
    numerator: DensePolynomial<F>,
    denominator: DensePolynomial<F>,
}

impl<F: Field> RationalFunction<F> {
    /// Builds a rational function from a numerator and denominator.
    pub fn new(
        numerator: DensePolynomial<F>,
        denominator: DensePolynomial<F>,
    ) -> Result<Self, FieldError> {
        let (numerator, denominator) = Self::canonicalize(numerator, denominator)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    /// Embeds a polynomial into `F(x)` as a denominator-one rational function.
    pub fn from_polynomial(polynomial: DensePolynomial<F>) -> Self {
        Self::new(polynomial, DensePolynomial::<F>::constant(F::one()))
            .expect("denominator one should define a valid rational function")
    }

    /// Builds the constant rational function with the given value.
    pub fn constant(value: F::Elem) -> Self {
        Self::from_polynomial(DensePolynomial::<F>::constant(value))
    }

    /// Returns the distinguished indeterminate `x`.
    pub fn indeterminate() -> Self {
        Self::from_polynomial(DensePolynomial::<F>::new(vec![F::zero(), F::one()]))
    }

    /// Returns the canonical numerator.
    pub fn numerator(&self) -> &DensePolynomial<F> {
        &self.numerator
    }

    /// Returns the canonical monic denominator.
    pub fn denominator(&self) -> &DensePolynomial<F> {
        &self.denominator
    }

    /// Returns whether the rational function is zero.
    pub fn is_zero(&self) -> bool {
        self.numerator.is_zero()
    }

    /// Returns whether the rational function equals `1`.
    pub fn is_one(&self) -> bool {
        self.numerator == self.denominator
    }

    /// Returns whether the rational function is constant in `F(x)`.
    ///
    /// Since the stored presentation is canonical with monic denominator, this
    /// is equivalent to asking whether both numerator and denominator have
    /// degree `0` (with the zero function treated as constant too).
    pub fn is_constant(&self) -> bool {
        self.denominator.degree() == Some(0) && self.numerator.degree().unwrap_or(0) == 0
    }

    /// Adds two rational functions and returns the canonical result.
    pub fn add(&self, rhs: &Self) -> Self {
        let numerator = self
            .numerator
            .mul(&rhs.denominator)
            .add(&rhs.numerator.mul(&self.denominator));
        let denominator = self.denominator.mul(&rhs.denominator);

        Self::new(numerator, denominator)
            .expect("product of non-zero denominators should stay non-zero")
    }

    /// Negates the rational function.
    pub fn neg(&self) -> Self {
        Self {
            numerator: self.numerator.neg(),
            denominator: self.denominator.clone(),
        }
    }

    /// Subtracts two rational functions and returns the canonical result.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// Multiplies two rational functions and returns the canonical result.
    pub fn mul(&self, rhs: &Self) -> Self {
        let numerator = self.numerator.mul(&rhs.numerator);
        let denominator = self.denominator.mul(&rhs.denominator);

        Self::new(numerator, denominator)
            .expect("product of non-zero denominators should stay non-zero")
    }

    /// Computes the multiplicative inverse when the rational function is
    /// non-zero.
    pub fn inverse(&self) -> Result<Self, FieldError> {
        if self.is_zero() {
            return Err(FieldError::DivisionByZero);
        }

        Self::new(self.denominator.clone(), self.numerator.clone())
    }

    /// Divides by another rational function when the divisor is non-zero.
    pub fn div(&self, rhs: &Self) -> Result<Self, FieldError> {
        Ok(self.mul(&rhs.inverse()?))
    }

    /// Returns the formal derivative of the rational function.
    ///
    /// `(p / q)' = (p' q - p q') / q^2`.
    pub fn derivative(&self) -> Self {
        let numerator = self
            .numerator
            .derivative()
            .mul(&self.denominator)
            .sub(&self.numerator.mul(&self.denominator.derivative()));
        let denominator = self.denominator.mul(&self.denominator);

        Self::new(numerator, denominator)
            .expect("square of a non-zero denominator should stay non-zero")
    }

    /// Raises the rational function to a nonnegative integer power.
    ///
    /// This uses binary exponentiation, so it performs `Θ(log exponent)`
    /// rational-function multiplications.
    pub fn pow_u64(&self, exponent: u64) -> Self {
        let mut result = Self::constant(F::one());
        let mut base = self.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            exp >>= 1;
            if exp > 0 {
                base = base.mul(&base);
            }
        }

        result
    }

    /// Raises the rational function to a nonnegative `u128` power.
    ///
    /// This stays crate-private because the current public educational surface
    /// only needs a smaller integer exponent API, while some Frobenius-side
    /// plumbing needs to work with field-order powers represented in `u128`.
    ///
    /// Like [`Self::pow_u64`], this uses binary exponentiation and therefore
    /// performs `Θ(log exponent)` rational-function multiplications.
    pub(crate) fn pow_u128(&self, exponent: u128) -> Self {
        let mut result = Self::constant(F::one());
        let mut base = self.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            exp >>= 1;
            if exp > 0 {
                base = base.mul(&base);
            }
        }

        result
    }

    fn canonicalize(
        numerator: DensePolynomial<F>,
        denominator: DensePolynomial<F>,
    ) -> Result<(DensePolynomial<F>, DensePolynomial<F>), FieldError> {
        if denominator.is_zero() {
            return Err(FieldError::DivisionByZero);
        }

        if numerator.is_zero() {
            return Ok((
                DensePolynomial::<F>::new(Vec::new()),
                DensePolynomial::<F>::constant(F::one()),
            ));
        }

        let gcd = numerator.gcd(&denominator);
        let numerator = numerator.quo(&gcd).map_err(Self::map_polynomial_error)?;
        let denominator = denominator.quo(&gcd).map_err(Self::map_polynomial_error)?;

        let denominator_leading = denominator
            .leading_coefficient()
            .expect("non-zero denominator has a leading coefficient");
        let denominator_scale = F::inverse(denominator_leading)?;

        Ok((
            numerator.scale(&denominator_scale),
            denominator.scale(&denominator_scale),
        ))
    }

    fn same_polynomials(lhs: &DensePolynomial<F>, rhs: &DensePolynomial<F>) -> bool {
        lhs == rhs
    }

    fn map_polynomial_error(error: PolynomialError) -> FieldError {
        match error {
            PolynomialError::DivisionByZeroPolynomial => FieldError::DivisionByZero,
            PolynomialError::NonInvertibleLeadingCoefficient => FieldError::NonInvertibleElement,
            _ => FieldError::Unsupported(
                "unexpected polynomial-domain error during rational-function normalization",
            ),
        }
    }
}

impl<F: Field> Clone for RationalFunction<F> {
    fn clone(&self) -> Self {
        Self {
            numerator: self.numerator.clone(),
            denominator: self.denominator.clone(),
        }
    }
}

impl<F: Field> fmt::Debug for RationalFunction<F>
where
    F::Elem: fmt::Debug,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RationalFunction")
            .field("numerator", &self.numerator.coefficients())
            .field("denominator", &self.denominator.coefficients())
            .finish()
    }
}

impl<F: Field> PartialEq for RationalFunction<F> {
    fn eq(&self, other: &Self) -> bool {
        Self::same_polynomials(&self.numerator, &other.numerator)
            && Self::same_polynomials(&self.denominator, &other.denominator)
    }
}

impl<F: FiniteField> PthRootExtraction for RationalFunction<F>
where
    DensePolynomial<F>: PthRootExtraction,
{
    /// Returns one `p`-th root of the rational function when it exists in `F(x)`.
    ///
    /// Because [`RationalFunction`] stores a reduced canonical presentation
    /// `numerator / denominator`, this implementation asks for `p`-th roots of
    /// those two coprime polynomials separately. Over a perfect field of
    /// characteristic `p`, that is exactly the condition for the rational
    /// function itself to be a `p`-th power in `F(x)`.
    fn pth_root(&self) -> Option<Self> {
        let numerator = self.numerator.pth_root()?;
        let denominator = self.denominator.pth_root()?;
        Self::new(numerator, denominator).ok()
    }
}

impl<F: FiniteField> RationalFunction<F> {
    /// Inverts the coordinate substitution `x' ↦ x^p` on a rational function
    /// when possible.
    ///
    /// This stays crate-private because it is not part of the generic public
    /// story of `F(x)`: it is a specialized helper for current Frobenius-side
    /// function-field plumbing.
    ///
    /// This helper is designed for inverting the current function-field
    /// pullback of absolute Frobenius `E -> E^(p)`. It is intentionally
    /// different from [`PthRootExtraction`]:
    ///
    /// - the coefficients stay fixed,
    /// - only the exponents are divided by `p`,
    /// - and the function must already lie in the image of substitution
    ///   `R(x') ↦ R(x^p)`.
    pub(crate) fn inverse_absolute_frobenius_pullback_from_twist(&self) -> Option<Self> {
        let numerator = self
            .numerator
            .inverse_absolute_frobenius_pullback_from_twist()?;
        let denominator = self
            .denominator
            .inverse_absolute_frobenius_pullback_from_twist()?;

        Self::new(numerator, denominator).ok()
    }
}
