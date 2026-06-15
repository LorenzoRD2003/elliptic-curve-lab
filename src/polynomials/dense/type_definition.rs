use crate::fields::traits::Field;

/// Dense polynomial over a field `F`, with coefficients stored in ascending
/// degree order.
///
/// This type is intentionally specialized to coefficients that live in a
/// field. That is more restrictive than the most general algebraic story, but
/// it matches the current scope of the project:
///
/// - the repository is currently centered on field-oriented algebra
/// - the user explicitly does not want to support coefficient rings that are
///   not fields yet
/// - adding polynomial arithmetic now is more useful than preserving maximum
///   abstract generality too early
///
/// In other words, this module currently models polynomials in `F[x]`, not the
/// more general `R[x]` for arbitrary rings `R`.
///
/// The storage convention is:
///
/// - `coefficients[0]` is the constant term
/// - `coefficients[1]` is the coefficient of `x`
/// - `coefficients[2]` is the coefficient of `x^2`
/// - and so on
///
/// For example, the vector `[a0, a1, a2]` represents
///
/// `a0 + a1*x + a2*x^2`
///
/// Important educational note:
///
/// this implementation trims trailing zero coefficients so the stored
/// representation stays canonical among dense vectors. As a consequence:
///
/// - `[5, 0, 0]` is normalized to `[5]`
/// - `[0, 0]` is normalized to `[]`
/// - [`DensePolynomial::degree`] reports the degree of the normalized dense
///   representation
#[derive(Debug)]
pub struct DensePolynomial<F: Field> {
    pub(super) coefficients: Vec<F::Elem>,
}

impl<F: Field> Clone for DensePolynomial<F> {
    fn clone(&self) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
        }
    }
}

impl<F: Field> PartialEq for DensePolynomial<F> {
    fn eq(&self, other: &Self) -> bool {
        self.coefficients.len() == other.coefficients.len()
            && self
                .coefficients
                .iter()
                .zip(&other.coefficients)
                .all(|(lhs, rhs)| F::eq(lhs, rhs))
    }
}

impl<F: Field> DensePolynomial<F> {
    /// Builds a dense polynomial from raw coefficients in ascending degree
    /// order.
    ///
    /// The constructor preserves the vector exactly as provided. It does not
    /// attempt any algebraic simplification beyond trimming trailing zero
    /// coefficients.
    pub fn new(coefficients: Vec<F::Elem>) -> Self {
        Self {
            coefficients: Self::trim_trailing_zero_coefficients(coefficients),
        }
    }

    /// Returns the stored coefficients in ascending degree order.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }

    /// Builds the constant polynomial with the given value.
    pub fn constant(value: F::Elem) -> Self {
        Self::new(vec![value])
    }

    /// Returns the number of stored coefficients.
    ///
    /// Because the representation trims trailing zeros, this is also the length
    /// of the canonical dense storage currently used by the type.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.coefficients.len()
    }

    /// Returns the degree of the normalized dense representation.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }

    /// Returns the leading stored coefficient, if any.
    ///
    /// This is simply the last coefficient in the normalized vector, so it is
    /// guaranteed to be non-zero whenever the polynomial is non-empty.
    pub fn leading_coefficient(&self) -> Option<&F::Elem> {
        self.coefficients.last()
    }

    /// Returns whether the polynomial is monic.
    ///
    /// A non-zero polynomial is monic when its leading coefficient equals the
    /// multiplicative identity of the coefficient field.
    ///
    /// The zero polynomial is intentionally not considered monic.
    pub fn is_monic(&self) -> bool {
        self.leading_coefficient()
            .is_some_and(|leading| F::eq(leading, &F::one()))
    }

    /// Returns the constant term, if any.
    ///
    /// This is the coefficient of `x^0`. For the canonical empty-vector
    /// representation of the zero polynomial, this returns `None`.
    pub fn constant_term(&self) -> Option<&F::Elem> {
        self.coefficients.first()
    }

    /// Returns whether the polynomial stores no coefficients at all.
    ///
    /// In this scaffold, an empty vector is allowed and is interpreted as a
    /// zero-polynomial representation with no explicit stored terms.
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// Removes trailing zero coefficients from a dense coefficient vector.
    ///
    /// This helper keeps the zero polynomial in the canonical empty-vector form
    /// and is used by the constructor so the representation remains normalized
    /// after creation and arithmetic.
    pub(super) fn trim_trailing_zero_coefficients(mut coefficients: Vec<F::Elem>) -> Vec<F::Elem> {
        while coefficients.last().is_some_and(F::is_zero) {
            coefficients.pop();
        }

        coefficients
    }

    /// Shifts a dense polynomial by multiplying it by `x^degree`.
    ///
    /// This helper is primarily used internally by division algorithms. The
    /// zero polynomial stays in its canonical empty-vector representation.
    pub(super) fn shift_by(&self, degree: usize) -> Self {
        if self.is_zero() {
            return Self::new(Vec::new());
        }

        let mut coefficients = vec![F::zero(); degree];
        coefficients.extend(self.coefficients.iter().cloned());
        Self::new(coefficients)
    }
}
