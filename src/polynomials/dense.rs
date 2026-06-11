use crate::fields::{Field, FiniteField, PthRootExtraction};
use crate::polynomials::{PolynomialError, SparsePolynomial, UnivariatePolynomial};

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
    coefficients: Vec<F::Elem>,
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
    fn trim_trailing_zero_coefficients(mut coefficients: Vec<F::Elem>) -> Vec<F::Elem> {
        while coefficients.last().is_some_and(F::is_zero) {
            coefficients.pop();
        }

        coefficients
    }

    /// Shifts a dense polynomial by multiplying it by `x^degree`.
    ///
    /// This helper is primarily used internally by division algorithms. The
    /// zero polynomial stays in its canonical empty-vector representation.
    fn shift_by(&self, degree: usize) -> Self {
        if self.is_zero() {
            return Self::new(Vec::new());
        }

        let mut coefficients = vec![F::zero(); degree];
        coefficients.extend(self.coefficients.iter().cloned());
        Self::new(coefficients)
    }

    /// Adds two dense polynomials coefficient-wise.
    ///
    /// The result length is the maximum of the two input lengths. Coefficients
    /// beyond the shorter input are treated as zero. The result is normalized
    /// through the constructor, so trailing zeros are trimmed away.
    pub fn add(&self, rhs: &Self) -> Self {
        let target_len = self.len().max(rhs.len());
        let mut coefficients = Vec::with_capacity(target_len);

        for index in 0..target_len {
            let lhs_coeff = self
                .coefficients
                .get(index)
                .cloned()
                .unwrap_or_else(F::zero);
            let rhs_coeff = rhs.coefficients.get(index).cloned().unwrap_or_else(F::zero);
            coefficients.push(F::add(&lhs_coeff, &rhs_coeff));
        }

        Self::new(coefficients)
    }

    /// Negates every coefficient of the polynomial.
    pub fn neg(&self) -> Self {
        Self::new(self.coefficients.iter().map(F::neg).collect())
    }

    /// Subtracts two dense polynomials coefficient-wise.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// Multiplies every coefficient of the polynomial by the same field
    /// element.
    pub fn scale(&self, scalar: &F::Elem) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .map(|coefficient| F::mul(coefficient, scalar))
                .collect(),
        )
    }

    /// Returns the formal derivative of the polynomial. If
    /// `f(x) = a0 + a1*x + a2*x^2 + ... + an*x^n`, then this
    ///  method returns `f'(x) = a1 + 2*a2*x + ... + n*an*x^(n-1)`.
    pub fn derivative(&self) -> Self {
        if self.len() <= 1 {
            return Self::new(Vec::new());
        }

        let coefficients = self
            .coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(degree, coefficient)| {
                let scalar = F::elem_from_u64(
                    u64::try_from(degree).expect("dense polynomial degree index should fit in u64"),
                );
                F::mul(coefficient, &scalar)
            })
            .collect();

        Self::new(coefficients)
    }

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

    /// Multiplies two dense polynomials using the naive quadratic algorithm.
    ///
    /// This implementation is intentionally straightforward and educational,
    /// and is not optimized.
    ///
    /// If either side is represented by an empty coefficient vector, the
    /// result keeps that empty zero representation.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::new(Vec::new());
        }

        let mut coefficients = vec![F::zero(); self.len() + rhs.len() - 1];

        for (lhs_degree, lhs_coeff) in self.coefficients.iter().enumerate() {
            for (rhs_degree, rhs_coeff) in rhs.coefficients.iter().enumerate() {
                let index = lhs_degree + rhs_degree;
                let term = F::mul(lhs_coeff, rhs_coeff);
                let updated = F::add(&coefficients[index], &term);
                coefficients[index] = updated;
            }
        }

        Self::new(coefficients)
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

impl<F: Field> UnivariatePolynomial<F> for DensePolynomial<F> {
    fn constant(value: F::Elem) -> Self {
        DensePolynomial::constant(value)
    }

    fn degree(&self) -> Option<usize> {
        DensePolynomial::degree(self)
    }

    fn leading_coefficient(&self) -> Option<&F::Elem> {
        DensePolynomial::leading_coefficient(self)
    }

    fn constant_term(&self) -> Option<&F::Elem> {
        DensePolynomial::constant_term(self)
    }

    fn is_zero(&self) -> bool {
        DensePolynomial::is_zero(self)
    }

    fn is_monic(&self) -> bool {
        DensePolynomial::is_monic(self)
    }

    fn add(&self, rhs: &Self) -> Self {
        DensePolynomial::add(self, rhs)
    }

    fn neg(&self) -> Self {
        DensePolynomial::neg(self)
    }

    fn sub(&self, rhs: &Self) -> Self {
        DensePolynomial::sub(self, rhs)
    }

    fn scale(&self, scalar: &F::Elem) -> Self {
        DensePolynomial::scale(self, scalar)
    }

    fn mul(&self, rhs: &Self) -> Self {
        DensePolynomial::mul(self, rhs)
    }

    fn derivative(&self) -> Self {
        DensePolynomial::derivative(self)
    }

    fn gcd(&self, rhs: &Self) -> Self {
        DensePolynomial::gcd(self, rhs)
    }

    fn make_monic(&self) -> Result<Self, PolynomialError> {
        DensePolynomial::make_monic(self)
    }
}

impl<F: FiniteField> PthRootExtraction for DensePolynomial<F>
where
    F::Elem: PthRootExtraction,
{
    fn pth_root(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::new(Vec::new()));
        }

        let characteristic = F::characteristic();
        let mut coefficients = Vec::new();

        for (degree, coefficient) in self.coefficients.iter().enumerate() {
            if F::is_zero(coefficient) {
                continue;
            }

            let degree_u64 = u64::try_from(degree).ok()?;
            if degree_u64 % characteristic != 0 {
                return None;
            }

            let root_degree = usize::try_from(degree_u64 / characteristic).ok()?;
            if coefficients.len() <= root_degree {
                coefficients.resize_with(root_degree + 1, F::zero);
            }

            coefficients[root_degree] = coefficient.pth_root()?;
        }

        Some(Self::new(coefficients))
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::fields::{Field, Fp, PthRootExtraction, Q};
    use crate::polynomials::{
        PolynomialError, SparsePolynomial, SparsePolynomialTerm, UnivariatePolynomial,
    };
    use crate::proptest_support::config::PolynomialStrategyConfig;
    use crate::proptest_support::polynomials::arb_dense_polynomial;

    use crate::polynomials::DensePolynomial;

    type F17 = Fp<17>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F17Sqrt3DensePthRootSpec,
        field: F17Sqrt3DensePthRoot,
        base: F17,
        non_residue: 3,
        name: "F17(sqrt(3)) for dense polynomial p-th-root tests",
    );

    fn f17_coefficients(values: &[u64]) -> Vec<<F17 as Field>::Elem> {
        values.iter().copied().map(F17::elem_from_u64).collect()
    }

    fn q_coefficients(values: &[(i64, i64)]) -> Vec<<Q as Field>::Elem> {
        values
            .iter()
            .map(|&(numerator, denominator)| {
                let numerator = Q::from_i64(numerator);
                let denominator = Q::from_i64(denominator);
                Q::div(&numerator, &denominator).expect("denominator should be non-zero")
            })
            .collect()
    }

    fn f17_sparse_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::elem_from_u64(coefficient),
            degree,
        }
    }

    #[test]
    fn dense_polynomial_preserves_storage_order_after_normalization() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[3, 15, 0, 7]));

        let coefficients = polynomial.coefficients();
        assert_eq!(coefficients.len(), 4);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(15)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(0)));
        assert!(F17::eq(&coefficients[3], &F17::elem_from_u64(7)));
        assert_eq!(polynomial.len(), 4);
        assert_eq!(polynomial.degree(), Some(3));
        assert!(F17::eq(
            polynomial
                .leading_coefficient()
                .expect("leading coefficient"),
            &F17::elem_from_u64(7)
        ));
    }

    #[test]
    fn dense_polynomial_allows_empty_storage_for_zero_representation() {
        let polynomial = DensePolynomial::<F17>::new(Vec::new());

        assert!(polynomial.is_zero());
        assert_eq!(polynomial.len(), 0);
        assert_eq!(polynomial.degree(), None);
        assert_eq!(polynomial.leading_coefficient(), None);
        assert_eq!(polynomial.constant_term(), None);
    }

    #[test]
    fn dense_polynomial_trims_trailing_zero_coefficients() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[5, 0, 0]));

        assert_eq!(polynomial.coefficients().len(), 1);
        assert_eq!(polynomial.degree(), Some(0));
        assert!(F17::eq(
            polynomial
                .leading_coefficient()
                .expect("leading coefficient"),
            &F17::elem_from_u64(5)
        ));
    }

    #[test]
    fn dense_polynomial_normalizes_all_zero_storage_to_empty() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[0, 0, 0]));

        assert!(polynomial.is_zero());
        assert_eq!(polynomial.coefficients(), &[]);
        assert_eq!(polynomial.degree(), None);
        assert_eq!(polynomial.leading_coefficient(), None);
    }

    #[test]
    fn dense_polynomial_single_coefficient_has_degree_zero() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[8]));

        assert_eq!(polynomial.degree(), Some(0));
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::elem_from_u64(8)
        ));
        assert!(F17::eq(
            polynomial
                .leading_coefficient()
                .expect("leading coefficient"),
            &F17::elem_from_u64(8)
        ));
        assert!(!polynomial.is_zero());
    }

    #[test]
    fn dense_polynomial_addition_is_coefficient_wise_over_f17() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[15, 14]));
        let sum = lhs.add(&rhs);

        let coefficients = sum.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(1)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(2)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(1)));
    }

    #[test]
    fn dense_polynomial_pth_root_over_prime_field_recovers_expected_coefficients() {
        let mut coefficients = vec![F17::zero(); 35];
        coefficients[0] = F17::elem_from_u64(4);
        coefficients[17] = F17::elem_from_u64(9);
        coefficients[34] = F17::elem_from_u64(3);
        let polynomial = DensePolynomial::<F17>::new(coefficients);

        let root = polynomial
            .pth_root()
            .expect("all non-zero term degrees are divisible by the characteristic");

        assert_eq!(
            root,
            DensePolynomial::<F17>::new(f17_coefficients(&[4, 9, 3]))
        );
        assert!(polynomial.has_pth_root());
    }

    #[test]
    fn dense_polynomial_pth_root_rejects_non_multiple_degree_terms() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 0, 0, 3]));

        assert_eq!(polynomial.pth_root(), None);
        assert!(!polynomial.has_pth_root());
    }

    #[test]
    fn dense_polynomial_pth_root_handles_the_zero_polynomial() {
        let polynomial = DensePolynomial::<F17>::new(Vec::new());

        assert_eq!(
            polynomial.pth_root(),
            Some(DensePolynomial::<F17>::new(Vec::new()))
        );
        assert!(polynomial.has_pth_root());
    }

    #[test]
    fn dense_polynomial_pth_root_uses_extension_field_coefficient_roots() {
        let generator = F17Sqrt3DensePthRoot::element(vec![F17::zero(), F17::one()]);
        let expected_root = DensePolynomial::<F17Sqrt3DensePthRoot>::new(vec![
            generator.clone(),
            F17Sqrt3DensePthRoot::one(),
        ]);

        let mut coefficients = vec![F17Sqrt3DensePthRoot::zero(); 18];
        coefficients[0] =
            F17Sqrt3DensePthRoot::pow(&generator, F17Sqrt3DensePthRoot::characteristic());
        coefficients[17] = F17Sqrt3DensePthRoot::one();
        let polynomial = DensePolynomial::<F17Sqrt3DensePthRoot>::new(coefficients);

        assert!(polynomial.pth_root() == Some(expected_root));
        assert!(
            polynomial
                .constant_term()
                .expect("constant term should be present")
                .pth_root()
                == Some(generator)
        );
    }

    #[test]
    fn dense_polynomial_multiplication_uses_naive_convolution_over_f17() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 4, 5]));
        let product = lhs.mul(&rhs);

        let coefficients = product.coefficients();
        assert_eq!(coefficients.len(), 4);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(10)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(13)));
        assert!(F17::eq(&coefficients[3], &F17::elem_from_u64(10)));
    }

    #[test]
    fn dense_polynomial_addition_works_over_q_too() {
        let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3)]));
        let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 3), (-2, 3), (5, 4)]));
        let sum = lhs.add(&rhs);

        let coefficients = sum.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(Q::eq(
            &coefficients[0],
            &Q::div(&Q::from_i64(5), &Q::from_i64(6)).unwrap()
        ));
        assert!(Q::eq(&coefficients[1], &Q::zero()));
        assert!(Q::eq(
            &coefficients[2],
            &Q::div(&Q::from_i64(5), &Q::from_i64(4)).unwrap()
        ));
    }

    #[test]
    fn dense_polynomial_constant_constructor_is_canonical() {
        let polynomial = DensePolynomial::<F17>::constant(F17::elem_from_u64(9));

        assert_eq!(polynomial.coefficients().len(), 1);
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::elem_from_u64(9)
        ));
    }

    #[test]
    fn dense_polynomial_negation_and_subtraction_work_over_f17() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[15, 14]));

        let neg_rhs = rhs.neg();
        assert_eq!(
            neg_rhs,
            DensePolynomial::<F17>::new(f17_coefficients(&[2, 3]))
        );

        let difference = lhs.sub(&rhs);
        assert_eq!(
            difference,
            DensePolynomial::<F17>::new(f17_coefficients(&[5, 8, 1]))
        );
    }

    #[test]
    fn dense_polynomial_scale_multiplies_every_coefficient() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5, 1]));
        let scaled = polynomial.scale(&F17::elem_from_u64(4));

        assert_eq!(
            scaled,
            DensePolynomial::<F17>::new(f17_coefficients(&[12, 3, 4]))
        );
    }

    #[test]
    fn dense_polynomial_derivative_drops_the_constant_term() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[4, 3, 5, 2]));

        assert_eq!(
            polynomial.derivative(),
            DensePolynomial::<F17>::new(f17_coefficients(&[3, 10, 6]))
        );
    }

    #[test]
    fn dense_polynomial_derivative_of_constant_is_zero() {
        let polynomial = DensePolynomial::<F17>::constant(F17::elem_from_u64(9));

        assert!(polynomial.derivative().is_zero());
        assert_eq!(polynomial.derivative().coefficients(), &[]);
    }

    #[test]
    fn dense_polynomial_derivative_trims_zero_tail_after_characteristic_cancellation() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 0, 0, 0, 0, 17]));

        assert!(polynomial.derivative().is_zero());
        assert_eq!(polynomial.derivative().coefficients(), &[]);
    }

    #[test]
    fn dense_polynomial_derivative_works_over_q_too() {
        let polynomial = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3), (3, 4)]));

        assert_eq!(
            polynomial.derivative(),
            DensePolynomial::<Q>::new(q_coefficients(&[(2, 3), (3, 2)]))
        );
    }

    #[test]
    fn dense_polynomial_manual_partial_eq_uses_field_equality() {
        let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(2, 4), (3, 6)]));
        let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (1, 2)]));
        let different = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (2, 3)]));

        assert_eq!(lhs, rhs);
        assert_ne!(lhs, different);
    }

    #[test]
    fn dense_polynomial_monic_helpers_work_over_f17() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[6, 3]));
        let monic = polynomial
            .make_monic()
            .expect("non-zero polynomial should normalize");

        assert!(!polynomial.is_monic());
        assert!(monic.is_monic());
        assert_eq!(
            monic,
            DensePolynomial::<F17>::new(f17_coefficients(&[2, 1]))
        );
    }

    #[test]
    fn dense_polynomial_zero_has_no_monic_normalization() {
        let polynomial = DensePolynomial::<F17>::new(Vec::new());
        let error = polynomial
            .make_monic()
            .expect_err("zero polynomial should not be monic-normalizable");

        assert_eq!(
            error,
            PolynomialError::ZeroPolynomialHasNoMonicNormalization
        );
        assert!(!polynomial.is_monic());
    }

    #[test]
    fn dense_polynomial_division_returns_zero_quotient_when_divisor_has_higher_degree() {
        let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5]));
        let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));

        let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should work");

        assert!(quotient.is_zero());
        assert_eq!(remainder, dividend);
    }

    #[test]
    fn dense_polynomial_division_handles_exact_division() {
        let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
        let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]));

        let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should work");

        assert_eq!(
            quotient,
            DensePolynomial::<F17>::new(f17_coefficients(&[2, 1]))
        );
        assert!(remainder.is_zero());
        assert_eq!(dividend, divisor.mul(&quotient).add(&remainder));
    }

    #[test]
    fn dense_polynomial_division_handles_non_zero_remainder() {
        let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 0, 1]));
        let divisor = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));

        let (quotient, remainder) = dividend.div_rem(&divisor).expect("division should work");

        assert_eq!(
            quotient,
            DensePolynomial::<F17>::new(f17_coefficients(&[0, 1]))
        );
        assert_eq!(
            remainder,
            DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]))
        );
        assert!(
            remainder.degree().expect("remainder should be non-zero") < divisor.degree().unwrap()
        );
        assert_eq!(dividend, divisor.mul(&quotient).add(&remainder));
    }

    #[test]
    fn dense_polynomial_division_by_constant_scales_coefficients() {
        let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[4, 8, 12]));
        let divisor = DensePolynomial::<F17>::constant(F17::elem_from_u64(4));

        let quotient = dividend.quo(&divisor).expect("division should work");
        let remainder = dividend.rem(&divisor).expect("division should work");

        assert_eq!(
            quotient,
            DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]))
        );
        assert!(remainder.is_zero());
    }

    #[test]
    fn dense_polynomial_division_rejects_zero_divisor() {
        let dividend = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]));
        let divisor = DensePolynomial::<F17>::new(Vec::new());

        let error = dividend
            .div_rem(&divisor)
            .expect_err("zero divisor should fail");
        assert_eq!(error, PolynomialError::DivisionByZeroPolynomial);
    }

    #[test]
    fn dense_polynomial_gcd_returns_monic_common_divisor() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 3, 3, 1]));

        let gcd = lhs.gcd(&rhs);

        assert_eq!(gcd, DensePolynomial::<F17>::new(f17_coefficients(&[1, 1])));
        assert!(gcd.is_monic());
    }

    #[test]
    fn dense_polynomial_gcd_of_coprimes_is_one() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 0, 1]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]));

        let gcd = lhs.gcd(&rhs);

        assert_eq!(gcd, DensePolynomial::<F17>::constant(F17::one()));
    }

    #[test]
    fn dense_polynomial_gcd_handles_zero_inputs() {
        let zero = DensePolynomial::<F17>::new(Vec::new());
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[2, 4]));

        assert_eq!(zero.gcd(&zero), DensePolynomial::<F17>::new(Vec::new()));
        assert_eq!(
            zero.gcd(&polynomial),
            DensePolynomial::<F17>::new(f17_coefficients(&[9, 1]))
        );
        assert_eq!(
            polynomial.gcd(&zero),
            DensePolynomial::<F17>::new(f17_coefficients(&[9, 1]))
        );
    }

    #[test]
    fn dense_polynomial_multiplication_works_over_q_too() {
        let lhs = DensePolynomial::<Q>::new(q_coefficients(&[(1, 2), (1, 3)]));
        let rhs = DensePolynomial::<Q>::new(q_coefficients(&[(2, 5), (3, 7)]));
        let product = lhs.mul(&rhs);

        let coefficients = product.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(Q::eq(
            &coefficients[0],
            &Q::div(&Q::from_i64(1), &Q::from_i64(5)).unwrap()
        ));
        assert!(Q::eq(
            &coefficients[1],
            &Q::div(&Q::from_i64(73), &Q::from_i64(210)).unwrap()
        ));
        assert!(Q::eq(
            &coefficients[2],
            &Q::div(&Q::from_i64(1), &Q::from_i64(7)).unwrap()
        ));
    }

    #[test]
    fn dense_polynomial_multiplication_preserves_empty_zero_representation() {
        let lhs = DensePolynomial::<F17>::new(Vec::new());
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2, 3]));
        let product = lhs.mul(&rhs);

        assert!(product.is_zero());
        assert_eq!(product.coefficients(), &[]);
    }

    #[test]
    fn dense_polynomial_addition_trims_zero_tail_after_cancellation() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 2]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[16, 15]));
        let sum = lhs.add(&rhs);

        assert!(sum.is_zero());
        assert_eq!(sum.coefficients(), &[]);
    }

    #[test]
    fn sparse_to_dense_conversion_fills_missing_degrees() {
        let sparse = SparsePolynomial::<F17>::new(vec![
            f17_sparse_term(3, 0),
            f17_sparse_term(5, 2),
            f17_sparse_term(1, 4),
        ]);
        let dense = DensePolynomial::<F17>::from(sparse);

        let coefficients = dense.coefficients();
        assert_eq!(coefficients.len(), 5);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[1], &F17::zero()));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(5)));
        assert!(F17::eq(&coefficients[3], &F17::zero()));
        assert!(F17::eq(&coefficients[4], &F17::elem_from_u64(1)));
    }

    #[test]
    fn sparse_to_dense_conversion_preserves_zero_polynomial() {
        let sparse = SparsePolynomial::<F17>::new(Vec::new());
        let dense = DensePolynomial::<F17>::from(sparse);

        assert!(dense.is_zero());
        assert_eq!(dense.coefficients(), &[]);
    }

    fn generic_addition<P>(lhs: &P, rhs: &P) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        lhs.add(rhs)
    }

    fn generic_derivative<P>(polynomial: &P) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        polynomial.derivative()
    }

    fn generic_gcd<P>(lhs: &P, rhs: &P) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        lhs.gcd(rhs)
    }

    #[test]
    fn dense_polynomial_implements_univariate_trait() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[3, 5]));
        let rhs = DensePolynomial::<F17>::constant(F17::elem_from_u64(1));
        let sum = generic_addition(&lhs, &rhs);

        assert_eq!(sum, DensePolynomial::<F17>::new(f17_coefficients(&[4, 5])));
        assert!(DensePolynomial::<F17>::constant(F17::one()).is_monic());
    }

    #[test]
    fn dense_polynomial_trait_derivative_uses_shared_surface() {
        let polynomial = DensePolynomial::<F17>::new(f17_coefficients(&[6, 5, 4]));

        assert_eq!(
            generic_derivative(&polynomial),
            DensePolynomial::<F17>::new(f17_coefficients(&[5, 8]))
        );
    }

    #[test]
    fn dense_polynomial_trait_gcd_uses_shared_surface() {
        let lhs = DensePolynomial::<F17>::new(f17_coefficients(&[2, 3, 1]));
        let rhs = DensePolynomial::<F17>::new(f17_coefficients(&[1, 3, 3, 1]));

        assert_eq!(
            generic_gcd(&lhs, &rhs),
            DensePolynomial::<F17>::new(f17_coefficients(&[1, 1]))
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(40))]

        #[test]
        fn property_dense_division_reconstructs_the_dividend(
            dividend in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
                max_len: 6,
                ..PolynomialStrategyConfig::default()
            }),
            divisor in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
                max_len: 4,
                ..PolynomialStrategyConfig::default()
            }),
        ) {
            prop_assume!(!divisor.is_zero());

            let (quotient, remainder) = dividend.div_rem(&divisor).expect("non-zero divisor should divide");
            prop_assert_eq!(divisor.mul(&quotient).add(&remainder), dividend);
            prop_assert!(remainder.is_zero() || remainder.degree().expect("non-zero remainder has a degree") < divisor.degree().expect("non-zero divisor has a degree"));
        }

        #[test]
        fn property_dense_gcd_divides_both_inputs(
            left in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
                max_len: 5,
                ..PolynomialStrategyConfig::default()
            }),
            right in arb_dense_polynomial::<F17>(PolynomialStrategyConfig {
                max_len: 5,
                ..PolynomialStrategyConfig::default()
            }),
        ) {
            let gcd = left.gcd(&right);
            if gcd.is_zero() {
                prop_assert!(left.is_zero() && right.is_zero());
            } else {
                prop_assert!(gcd.is_monic());
                prop_assert!(left.rem(&gcd).expect("non-zero gcd should divide the left operand").is_zero());
                prop_assert!(right.rem(&gcd).expect("non-zero gcd should divide the right operand").is_zero());
            }
        }
    }
}
