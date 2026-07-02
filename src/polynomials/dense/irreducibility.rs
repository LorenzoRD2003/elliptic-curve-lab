use crate::polynomials::{
    DensePolynomial, PolynomialError,
    irreducibility::{IrreducibilityBackend, IrreducibilityStatus},
};

impl<F: IrreducibilityBackend> DensePolynomial<F> {
    /// Returns a structured irreducibility classification for this dense
    /// polynomial.
    ///
    /// The exact backend used depends on the coefficient field family:
    ///
    /// - prime fields `Fp<M, LIMBS>` currently use an exhaustive educational search
    ///   over monic candidate divisors
    /// - algebraically closed backends such as `ComplexApprox` can conclude
    ///   that every degree-`>= 2` polynomial is reducible, even when no
    ///   explicit factorization witness is currently returned
    /// - `Q` currently uses an exact but partial backend that either
    ///   certifies an answer or returns a typed inconclusive error
    ///
    /// Classification conventions:
    ///
    /// - degree `0` polynomials are reported as
    ///   [`IrreducibilityStatus::Constant`]
    /// - degree `1` polynomials are reported as
    ///   [`IrreducibilityStatus::Linear`]
    /// - constants are not considered irreducible
    /// - linears are considered irreducible
    ///
    /// TODO:
    /// - add a Rabin-style irreducibility test for finite fields once the
    ///   supporting polynomial infrastructure is mature enough
    /// - extend the exact partial backend for `Q` into a complete decision
    ///   procedure once integer-polynomial factorization infrastructure exists
    pub fn irreducibility_status(&self) -> Result<IrreducibilityStatus<F>, PolynomialError> {
        F::irreducibility_status_impl(self)
    }

    /// Returns whether this dense polynomial is irreducible in the current
    /// backend.
    ///
    /// This is the boolean convenience wrapper around
    /// [`DensePolynomial::irreducibility_status`].
    pub fn is_irreducible(&self) -> Result<bool, PolynomialError> {
        Ok(matches!(
            self.irreducibility_status()?,
            IrreducibilityStatus::Linear | IrreducibilityStatus::Irreducible
        ))
    }
}
