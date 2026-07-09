use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::fields::{Q, traits::Field};
use crate::numerics::gcd_bigint;
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
}
