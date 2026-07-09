/// Exponent vector of a multivariate monomial.
///
/// If the ambient variables are ordered as
///
/// `x_0, x_1, ..., x_{n-1}`
///
/// then the exponent vector `[e0, e1, ..., e_{n-1}]` represents
///
/// `x_0^e0 * x_1^e1 * ... * x_{n-1}^e_{n-1}`
///
/// For example, in arity `3`:
///
/// - `[2, 0, 1]` represents `x_0^2 * x_2`
/// - `[0, 1, 0]` represents `x_1`
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Monomial {
    exponents: Vec<usize>,
}

impl Monomial {
    /// Builds a monomial from an explicit exponent vector.
    pub fn new(exponents: Vec<usize>) -> Self {
        Self { exponents }
    }

    /// Returns the ambient number of variables.
    pub fn arity(&self) -> usize {
        self.exponents.len()
    }

    /// Returns the exponent vector in the ordered variable list.
    pub fn exponents(&self) -> &[usize] {
        &self.exponents
    }

    /// Returns the total degree of the monomial.
    pub fn total_degree(&self) -> usize {
        self.exponents.iter().sum()
    }

    /// Multiplies two monomials of the same arity by adding exponents
    /// component-wise.
    pub fn mul(&self, rhs: &Self) -> Option<Self> {
        if self.arity() != rhs.arity() {
            return None;
        }

        let exponents = self
            .exponents
            .iter()
            .zip(rhs.exponents())
            .map(|(lhs, rhs)| lhs + rhs)
            .collect();

        Some(Self { exponents })
    }
}
