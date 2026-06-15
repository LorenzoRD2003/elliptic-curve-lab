use crate::fields::traits::{Field, FiniteField, PthRootExtraction};
use crate::polynomials::{
    SparsePolynomial, sparse::SparsePolynomialTerm, traits::UnivariatePolynomial,
};

impl<F: Field> UnivariatePolynomial<F> for SparsePolynomial<F> {
    fn constant(value: F::Elem) -> Self {
        SparsePolynomial::constant(value)
    }

    fn degree(&self) -> Option<usize> {
        SparsePolynomial::degree(self)
    }

    fn leading_coefficient(&self) -> Option<&F::Elem> {
        SparsePolynomial::leading_coefficient(self)
    }

    fn constant_term(&self) -> Option<&F::Elem> {
        SparsePolynomial::constant_term(self)
    }

    fn is_zero(&self) -> bool {
        SparsePolynomial::is_zero(self)
    }

    fn add(&self, rhs: &Self) -> Self {
        SparsePolynomial::add(self, rhs)
    }

    fn neg(&self) -> Self {
        SparsePolynomial::neg(self)
    }

    fn sub(&self, rhs: &Self) -> Self {
        SparsePolynomial::sub(self, rhs)
    }

    fn scale(&self, scalar: &F::Elem) -> Self {
        SparsePolynomial::scale(self, scalar)
    }

    fn mul(&self, rhs: &Self) -> Self {
        SparsePolynomial::mul(self, rhs)
    }

    fn derivative(&self) -> Self {
        SparsePolynomial::derivative(self)
    }

    fn gcd(&self, rhs: &Self) -> Self {
        SparsePolynomial::gcd(self, rhs)
    }
}

impl<F: FiniteField> PthRootExtraction for SparsePolynomial<F>
where
    F::Elem: PthRootExtraction,
{
    fn pth_root(&self) -> Option<Self> {
        let characteristic = F::characteristic();
        let mut terms = Vec::with_capacity(self.terms.len());

        for term in &self.terms {
            let degree_u64 = u64::try_from(term.degree).ok()?;
            if degree_u64 % characteristic != 0 {
                return None;
            }

            terms.push(SparsePolynomialTerm {
                coefficient: term.coefficient.pth_root()?,
                degree: usize::try_from(degree_u64 / characteristic).ok()?,
            });
        }

        Some(Self::new(terms))
    }
}
