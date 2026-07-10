use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::BinaryQuadraticFormError;
use crate::fields::{Q, traits::Field};
use crate::numerics::{floor_div_bigint_by_positive, gcd_bigint};
use crate::polynomials::{MultivariatePolynomial, multivariate::MultivariateTerm};

/// Integral binary quadratic form `ax² + bxy + cy²`.
///
/// The stored data is the integral ternary `(a, b, c)`. The polynomial view is
/// derived by embedding those integers into `Q`, so polynomial operations can
/// be reused without weakening the public invariant that the form itself has
/// integral coefficients.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinaryQuadraticForm {
    a: BigInt,
    b: BigInt,
    c: BigInt,
}

impl BinaryQuadraticForm {
    /// Builds the integral form `ax² + bxy + cy²`.
    ///
    /// This constructor accepts every integral ternary. Mathematical
    /// restrictions used by the class-group story, such as primitiveness or
    /// positive-definiteness, are exposed as separate predicates.
    ///
    /// Complexity: `Θ(1)` ownership moves.
    pub fn new(a: BigInt, b: BigInt, c: BigInt) -> Self {
        Self { a, b, c }
    }

    /// Builds the integral form with prescribed leading coefficient, middle
    /// coefficient, and discriminant.
    ///
    /// This is the shared internal inverse to `Δ = b² − 4ac`: it computes
    /// `c = (b² − Δ)/(4a)` and returns `None` when `a = 0` or that quotient is
    /// not integral. Callers remain responsible for any semantic checks such
    /// as primitiveness, positive-definiteness, or reducedness.
    ///
    /// Complexity: exact big-integer multiplication, subtraction, one
    /// divisibility check, and one exact division.
    pub(crate) fn from_leading_middle_discriminant(
        leading: BigInt,
        middle: BigInt,
        discriminant: &BigInt,
    ) -> Option<Self> {
        let denominator = BigInt::from(4u8) * &leading;
        if denominator.is_zero() {
            return None;
        }

        let numerator = &middle * &middle - discriminant;

        if (&numerator % &denominator) != BigInt::zero() {
            return None;
        }

        Some(Self::new(leading, middle, numerator / denominator))
    }

    /// Returns the coefficient of `x²`.
    pub fn a(&self) -> &BigInt {
        &self.a
    }

    /// Returns the coefficient of `xy`.
    pub fn b(&self) -> &BigInt {
        &self.b
    }

    /// Returns the coefficient of `y²`.
    pub fn c(&self) -> &BigInt {
        &self.c
    }

    /// Returns the integral coefficient ternary `(a, b, c)`.
    pub fn coefficients(&self) -> (&BigInt, &BigInt, &BigInt) {
        (&self.a, &self.b, &self.c)
    }

    /// Builds the multivariate-polynomial view of this form over `Q`.
    ///
    /// The returned polynomial has arity `2`, with variables ordered as
    /// `x = x₀` and `y = x₁`. Its terms are exactly the nonzero terms among
    /// `ax²`, `bxy`, and `cy²`, normalized by [`MultivariatePolynomial`].
    ///
    /// Cost: three exact integer embeddings into `Q` and normalization of at
    /// most three candidate terms.
    pub fn polynomial(&self) -> MultivariatePolynomial<Q> {
        MultivariatePolynomial::<Q>::new(
            2,
            vec![
                MultivariateTerm::from_exponents(Q::from_bigint(&self.a), vec![2, 0]),
                MultivariateTerm::from_exponents(Q::from_bigint(&self.b), vec![1, 1]),
                MultivariateTerm::from_exponents(Q::from_bigint(&self.c), vec![0, 2]),
            ],
        )
        .expect("binary quadratic form monomials have arity 2")
    }

    /// Returns the discriminant `Δ = b² − 4ac`.
    ///
    /// Cost: exact big-integer arithmetic in the coefficient sizes.
    pub fn discriminant(&self) -> BigInt {
        (&self.b * &self.b) - BigInt::from(4u8) * &self.a * &self.c
    }

    /// Returns whether `gcd(a, b, c) = 1`.
    ///
    /// The zero ternary is not primitive because its gcd is `0`.
    ///
    /// Complexity: two exact integer gcd computations.
    pub fn is_primitive(&self) -> bool {
        let gcd_ab = gcd_bigint(&self.a, &self.b);
        gcd_bigint(&gcd_ab, &self.c).is_one()
    }

    /// Returns whether this form is positive definite.
    ///
    /// For an integral binary quadratic form this is equivalent to
    /// `a > 0` and `Δ < 0`.
    ///
    /// Complexity: one discriminant computation and two sign checks.
    pub fn is_positive_definite(&self) -> bool {
        self.a > BigInt::zero() && self.discriminant() < BigInt::zero()
    }

    /// Returns whether this positive-definite form is reduced.
    ///
    /// The convention used here is:
    /// - `|b| ≤ a ≤ c`
    /// - if `|b| = a` or `a = c`, then `b ≥ 0`.
    ///
    ///
    /// Forms that are not positive definite also return `false`.
    pub fn is_reduced_positive_definite(&self) -> bool {
        self.is_positive_definite() && self.satisfies_positive_definite_reduced_inequalities()
    }

    /// Reduces a positive-definite binary quadratic form.
    ///
    /// This is Gauss reduction for integral positive-definite forms. It
    /// preserves the discriminant and returns the reduced representative
    /// satisfying [`Self::is_reduced_positive_definite`].
    ///
    /// The algorithm repeatedly applies two unimodular changes of variables.
    /// First it normalizes the middle coefficient by substituting
    /// `x ↦ x + qy`, which sends
    ///
    /// (a, b, c) ↦ (a, b + 2aq, aq² + bq + c).
    ///
    /// The chosen integer `q` puts the new middle coefficient in the range
    /// `−a < b ≤ a`. If the right coefficient is then smaller than the left
    /// one, or if `a = c` with `b < 0`, it swaps the variables by
    /// `(x, y) ↦ (−y, x)`, which sends
    ///
    /// (a, b, c) ↦ (c, −b, a).
    ///
    /// Both transformations preserve `Δ = b² − 4ac`.
    pub fn reduce_positive_definite(&self) -> Result<Self, BinaryQuadraticFormError> {
        if !self.is_positive_definite() {
            return Err(BinaryQuadraticFormError::NotPositiveDefinite);
        }

        let mut reduced = self.clone();
        while !reduced.is_reduced_positive_definite() {
            reduced.normalize_middle_coefficient();

            if reduced.a > reduced.c || (reduced.a == reduced.c && reduced.b < BigInt::zero()) {
                reduced = Self::new(reduced.c, -reduced.b, reduced.a);
            }
        }

        Ok(reduced)
    }

    /// Returns the conjugate form `(a, −b, c)`.
    ///
    /// This is the form-side shadow of conjugating the corresponding ideal
    /// class in the later class-group story.
    ///
    /// Cost: exact big-integer clones and one negation.
    pub fn conjugate(&self) -> Self {
        Self::new(self.a.clone(), -&self.b, self.c.clone())
    }

    /// Evaluates `ax² + bxy + cy²` at integral coordinates.
    ///
    /// Cost: exact big-integer arithmetic in the coefficient and input sizes.
    pub fn evaluate_integral(&self, x: &BigInt, y: &BigInt) -> BigInt {
        &self.a * x * x + &self.b * x * y + &self.c * y * y
    }

    fn normalize_middle_coefficient(&mut self) {
        let two_a = BigInt::from(2u8) * &self.a;
        let q = floor_div_bigint_by_positive(&(&self.a - &self.b), &two_a);

        if q.is_zero() {
            return;
        }

        self.c = &self.a * &q * &q + &self.b * &q + &self.c;
        self.b = &self.b + two_a * q;
    }

    fn satisfies_positive_definite_reduced_inequalities(&self) -> bool {
        let abs_b = self.b.abs();
        abs_b <= self.a
            && self.a <= self.c
            && (abs_b != self.a || self.b >= BigInt::zero())
            && (self.a != self.c || self.b >= BigInt::zero())
    }
}
