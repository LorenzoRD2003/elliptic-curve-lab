use crate::fields::traits::Field;
use crate::polynomials::{DensePolynomial, SparsePolynomial, sparse::SparsePolynomialTerm};

impl<F: Field> SparsePolynomial<F> {
    /// Adds two sparse polynomials.
    ///
    /// The implementation concatenates the term lists and lets the constructor
    /// perform normalization via [`SparsePolynomial::normalize_terms`]. This
    /// is simple and sufficient for the current educational stage of the
    /// project.
    pub fn add(&self, rhs: &Self) -> Self {
        let mut terms = Vec::with_capacity(self.len() + rhs.len());
        terms.extend(self.terms.iter().cloned());
        terms.extend(rhs.terms.iter().cloned());
        Self::new(terms)
    }

    /// Negates every stored coefficient of the polynomial.
    pub fn neg(&self) -> Self {
        Self::new(
            self.terms
                .iter()
                .map(|term| SparsePolynomialTerm {
                    coefficient: F::neg(&term.coefficient),
                    degree: term.degree,
                })
                .collect(),
        )
    }

    /// Subtracts two sparse polynomials.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// Multiplies every stored coefficient by the same field element.
    pub fn scale(&self, scalar: &F::Elem) -> Self {
        Self::new(
            self.terms
                .iter()
                .map(|term| SparsePolynomialTerm {
                    coefficient: F::mul(&term.coefficient, scalar),
                    degree: term.degree,
                })
                .collect(),
        )
    }

    /// Multiplies two sparse polynomials using the naive pairwise product.
    ///
    /// Every term on the left is multiplied by every term on the right. The
    /// constructor is then used to normalize the resulting term list via
    /// [`SparsePolynomial::normalize_terms`].
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_empty() || rhs.is_empty() {
            return Self::new(Vec::new());
        }

        let mut terms = Vec::with_capacity(self.len() * rhs.len());

        for lhs_term in &self.terms {
            for rhs_term in &rhs.terms {
                terms.push(SparsePolynomialTerm {
                    coefficient: F::mul(&lhs_term.coefficient, &rhs_term.coefficient),
                    degree: lhs_term.degree + rhs_term.degree,
                });
            }
        }

        Self::new(terms)
    }

    /// Returns the formal derivative of the polynomial.
    ///
    /// Each term `a*x^d` with `d > 0` becomes `d*a*x^(d-1)`. Constant terms
    /// disappear, and any coefficient that becomes zero in the base field is
    /// removed by the sparse constructor's normalization.
    pub fn derivative(&self) -> Self {
        let terms = self
            .terms
            .iter()
            .filter(|term| term.degree > 0)
            .map(|term| SparsePolynomialTerm {
                coefficient: F::mul(
                    &term.coefficient,
                    &F::elem_from_u64(
                        u64::try_from(term.degree)
                            .expect("sparse polynomial degree should fit in u64"),
                    ),
                ),
                degree: term.degree - 1,
            })
            .collect();

        Self::new(terms)
    }

    /// Computes a greatest common divisor by delegating to the current dense
    /// Euclidean algorithm and converting the monic result back to sparse
    /// storage.
    ///
    /// This keeps the shared algebraic answer available for sparse
    /// polynomials without yet duplicating polynomial division machinery in
    /// the sparse representation.
    pub fn gcd(&self, rhs: &Self) -> Self {
        let left_dense = DensePolynomial::<F>::from(self.clone());
        let right_dense = DensePolynomial::<F>::from(rhs.clone());
        let gcd_dense = left_dense.gcd(&right_dense);
        Self::from(gcd_dense)
    }
}
