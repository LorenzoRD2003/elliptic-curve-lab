use crate::fields::traits::Field;
use crate::polynomials::{DensePolynomial, SparsePolynomial, sparse::SparsePolynomialTerm};

impl<F: Field> From<DensePolynomial<F>> for SparsePolynomial<F> {
    /// Converts a dense polynomial into its canonical sparse term list.
    ///
    /// Zero coefficients are discarded, so only the non-zero degrees survive
    /// in the sparse representation.
    fn from(polynomial: DensePolynomial<F>) -> Self {
        let terms = polynomial
            .coefficients()
            .iter()
            .enumerate()
            .filter_map(|(degree, coefficient)| {
                if F::is_zero(coefficient) {
                    None
                } else {
                    Some(SparsePolynomialTerm {
                        coefficient: coefficient.clone(),
                        degree,
                    })
                }
            })
            .collect();

        Self::new(terms)
    }
}
