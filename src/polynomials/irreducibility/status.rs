use crate::fields::Field;
use crate::polynomials::DensePolynomial;

/// Explanation for a reducibility result when the library does not currently
/// return an explicit factorization witness.
///
/// This keeps the educational API honest: some backends can certify
/// reducibility from field-theoretic metadata alone, even when the crate is
/// not yet set up to produce concrete factors in that backend.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReducibilityReason {
    /// The base field is algebraically closed, so every polynomial of degree
    /// at least `2` has a non-trivial factorization.
    AlgebraicallyClosedField,
}

/// Classification of a dense univariate polynomial over a field with respect
/// to irreducibility.
///
/// This enum is intentionally educational:
///
/// - `Constant` and `Linear` are called out explicitly instead of being folded
///   into `Irreducible` or `Reducible`
/// - `Reducible` carries a witness factorization when the current backend can
///   compute one
/// - `ReducibleWithoutWitness` allows the library to report a mathematically
///   correct reducibility conclusion without inventing unsupported factors
///
/// In the current implementation, the witness returned by
/// `irreducibility_status` satisfies
///
/// `polynomial = divisor * quotient`
///
/// in the original field, not merely after monic normalization.
#[derive(Clone, Debug)]
pub enum IrreducibilityStatus<F: Field> {
    /// The polynomial has degree `0` or is the zero polynomial.
    Constant,
    /// The polynomial has degree `1`.
    Linear,
    /// No non-trivial factor was found by the current backend algorithm.
    Irreducible,
    /// A non-trivial factorization was found.
    Reducible {
        /// A non-constant proper divisor of the polynomial.
        divisor: DensePolynomial<F>,
        /// The quotient satisfying `polynomial = divisor * quotient`.
        quotient: DensePolynomial<F>,
    },
    /// The polynomial is known to be reducible, but the current backend does
    /// not expose an explicit factorization witness.
    ReducibleWithoutWitness {
        /// High-level explanation for the reducibility conclusion.
        reason: ReducibilityReason,
    },
}

impl<F: Field> PartialEq for IrreducibilityStatus<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Constant, Self::Constant)
            | (Self::Linear, Self::Linear)
            | (Self::Irreducible, Self::Irreducible) => true,
            (
                Self::Reducible {
                    divisor: lhs_divisor,
                    quotient: lhs_quotient,
                },
                Self::Reducible {
                    divisor: rhs_divisor,
                    quotient: rhs_quotient,
                },
            ) => lhs_divisor == rhs_divisor && lhs_quotient == rhs_quotient,
            (
                Self::ReducibleWithoutWitness { reason: lhs_reason },
                Self::ReducibleWithoutWitness { reason: rhs_reason },
            ) => lhs_reason == rhs_reason,
            _ => false,
        }
    }
}
