use crate::fields::traits::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};

impl<F: Field> DensePolynomial<F> {
    /// Returns the monic normalization of the polynomial.
    ///
    /// For a non-zero polynomial `f(x)`, this divides every coefficient by the
    /// leading coefficient so the result has leading coefficient `1`.
    ///
    /// The zero polynomial has no monic normalization and is therefore
    /// rejected.
    pub fn make_monic(&self) -> Result<Self, PolynomialError> {
        let Some(leading) = self.leading_coefficient() else {
            return Err(PolynomialError::ZeroPolynomialHasNoMonicNormalization);
        };

        let inverse =
            F::inverse(leading).map_err(|_| PolynomialError::NonInvertibleLeadingCoefficient)?;
        Ok(self.scale(&inverse))
    }

    /// Divides one dense polynomial by another and returns the quotient and
    /// remainder.
    ///
    /// This method implements the classical univariate Euclidean division
    /// algorithm over a field. It returns polynomials `q(x)` and `r(x)` such
    /// that `self = divisor * q + r`
    /// with either `r = 0`, or `deg(r) < deg(divisor)`.
    ///
    /// The divisor must be non-zero.
    pub fn div_rem(&self, divisor: &Self) -> Result<(Self, Self), PolynomialError> {
        if divisor.is_zero() {
            return Err(PolynomialError::DivisionByZeroPolynomial);
        }

        if self.is_zero() {
            return Ok((Self::new(Vec::new()), Self::new(Vec::new())));
        }

        let Some(divisor_degree) = divisor.degree() else {
            return Err(PolynomialError::DivisionByZeroPolynomial);
        };
        let divisor_leading = divisor
            .leading_coefficient()
            .expect("non-zero divisor has a leading coefficient")
            .clone();

        if self.degree().expect("non-zero dividend has a degree") < divisor_degree {
            return Ok((Self::new(Vec::new()), self.clone()));
        }

        let mut quotient = Self::new(Vec::new());
        let mut remainder = self.clone();

        while let Some(remainder_degree) = remainder.degree() {
            if remainder_degree < divisor_degree {
                break;
            }

            let remainder_leading = remainder
                .leading_coefficient()
                .expect("non-zero remainder has a leading coefficient")
                .clone();
            let degree_gap = remainder_degree - divisor_degree;
            let scale = F::div(&remainder_leading, &divisor_leading)
                .map_err(|_| PolynomialError::NonInvertibleLeadingCoefficient)?;

            let quotient_term = Self::constant(scale.clone()).shift_by(degree_gap);
            quotient = quotient.add(&quotient_term);

            let subtraction_term = divisor.scale(&scale).shift_by(degree_gap);
            remainder = remainder.sub(&subtraction_term);
        }

        Ok((quotient, remainder))
    }

    /// Returns only the quotient of Euclidean division.
    pub fn quo(&self, divisor: &Self) -> Result<Self, PolynomialError> {
        let (quotient, _) = self.div_rem(divisor)?;
        Ok(quotient)
    }

    /// Returns only the remainder of Euclidean division.
    pub fn rem(&self, divisor: &Self) -> Result<Self, PolynomialError> {
        let (_, remainder) = self.div_rem(divisor)?;
        Ok(remainder)
    }

    /// Computes the greatest common divisor of two univariate polynomials over
    /// a field and returns it in monic form whenever it is non-zero.
    ///
    /// The implementation uses the classical Euclidean algorithm on
    /// remainders:
    ///
    /// `gcd(a, b) = gcd(b, a mod b)`
    ///
    /// and stops when the remainder becomes zero.
    ///
    /// Normalization policy:
    ///
    /// - if at least one input is non-zero, the returned gcd is monic
    /// - if both inputs are zero, the returned gcd is the zero polynomial
    ///
    /// This choice keeps the result canonical over a field, where gcds are
    /// only defined up to multiplication by a non-zero scalar.
    pub fn gcd(&self, rhs: &Self) -> Self {
        if self.is_zero() && rhs.is_zero() {
            return Self::new(Vec::new());
        }

        if self.is_zero() {
            return rhs
                .make_monic()
                .expect("non-zero polynomial admits monic normalization");
        }

        if rhs.is_zero() {
            return self
                .make_monic()
                .expect("non-zero polynomial admits monic normalization");
        }

        let mut a = self.clone();
        let mut b = rhs.clone();

        while !b.is_zero() {
            let remainder = a
                .rem(&b)
                .expect("euclidean step divides by a non-zero polynomial");
            a = b;
            b = remainder;
        }

        a.make_monic()
            .expect("final non-zero gcd admits monic normalization")
    }
}
