use crate::fields::traits::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};

impl<F: Field> DensePolynomial<F> {
    /// Interpolates a univariate polynomial over a field from sample points using
    /// the classical Lagrange formula.
    ///
    /// Given distinct sample points `(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})`
    /// this method constructs the unique polynomial `p(x)` of degree at most
    /// `n - 1` such that `p(x_i) = y_i` for every sample index `i`.
    ///
    /// The implementation is intentionally direct and educational:
    ///
    /// - it builds each Lagrange basis polynomial explicitly
    /// - it uses the already available dense polynomial arithmetic
    /// - it favors clarity over asymptotic efficiency
    ///
    /// Behavior notes:
    ///
    /// - an empty sample list returns the zero polynomial
    /// - a single sample returns the corresponding constant polynomial
    /// - repeated `x` coordinates are rejected
    ///
    /// TODO:
    ///
    /// - add Newton interpolation for an incremental formulation
    /// - add barycentric interpolation for numerically lighter repeated evaluation
    /// - add subproduct-tree / divide-and-conquer interpolation once performance
    ///   becomes relevant
    pub fn lagrange_interpolate(samples: &[(F::Elem, F::Elem)]) -> Result<Self, PolynomialError> {
        if samples.is_empty() {
            return Ok(Self::new(Vec::new()));
        }

        let mut result = Self::new(Vec::new());

        for (i, (x_i, y_i)) in samples.iter().enumerate() {
            let mut numerator = Self::new(vec![F::one()]);
            let mut denominator = F::one();

            for (j, (x_j, _)) in samples.iter().enumerate() {
                if i == j {
                    continue;
                }

                if F::eq(x_i, x_j) {
                    return Err(PolynomialError::DuplicateInterpolationAbscissa);
                }

                let factor = DensePolynomial::<F>::new(vec![F::neg(x_j), F::one()]);
                numerator = numerator.mul(&factor);

                let difference = F::sub(x_i, x_j);
                denominator = F::mul(&denominator, &difference);
            }

            let scaling = F::div(y_i, &denominator)
                .map_err(|_| PolynomialError::NonInvertibleInterpolationDenominator)?;
            let scaled_basis = numerator.scale(&scaling);
            result = result.add(&scaled_basis);
        }

        Ok(result)
    }
}
