use crate::fields::traits::*;
use crate::polynomials::{DensePolynomial, SparsePolynomial};

impl<F: Field> From<SparsePolynomial<F>> for DensePolynomial<F> {
    /// Converts a sparse polynomial into its canonical dense coefficient
    /// vector.
    ///
    /// Missing degrees are filled with zero coefficients, and the final dense
    /// representation is normalized through [`DensePolynomial::new`].
    fn from(polynomial: SparsePolynomial<F>) -> Self {
        let Some(max_degree) = polynomial.degree() else {
            return Self::new(Vec::new());
        };

        let mut coefficients = vec![F::zero(); max_degree + 1];

        for term in polynomial.terms() {
            coefficients[term.degree] = term.coefficient.clone();
        }

        Self::new(coefficients)
    }
}
