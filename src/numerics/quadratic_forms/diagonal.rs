use num_bigint::BigUint;
use num_traits::Zero;

use crate::numerics::{
    cornacchia::cornacchia_primitive_solutions,
    quadratic_forms::{DiagonalBinaryQuadraticRepresentation, QuadraticFormError},
};

/// The positive diagonal binary quadratic form `x² + d y²`.
///
/// This type is a conceptual layer over the current Cornacchia-backed
/// implementation. It lets callers ask representation questions about a form
/// directly, without depending on the algorithm used internally.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiagonalBinaryQuadraticForm {
    d: BigUint,
}

impl DiagonalBinaryQuadraticForm {
    /// Builds the diagonal form `x² + d y²`.
    ///
    /// The current type models positive diagonal forms, so `d = 0` is rejected.
    ///
    /// Complexity: `Θ(1)`.
    pub fn new(d: BigUint) -> Result<Self, QuadraticFormError> {
        if d.is_zero() {
            return Err(QuadraticFormError::ZeroDiagonalCoefficient);
        }
        Ok(Self { d })
    }

    /// Returns the coefficient `d` in `x² + d y²`.
    pub fn coefficient(&self) -> &BigUint {
        &self.d
    }

    /// Returns the primitive non-negative representations of `m` by this form.
    ///
    /// A representation is primitive when `gcd(x, y) = 1`. Internally this
    /// delegates to Cornacchia over all modular square roots of `-d` modulo
    /// `m`, then maps the resulting algorithm-level solutions into this
    /// form-level value object.
    ///
    /// For non-square-free `m`, this method intentionally exposes only the
    /// primitive representation surface; non-primitive representation
    /// enumeration is left to a later, explicit API.
    ///
    /// Complexity: dominated by the modular square-root computation for
    /// `-d mod m` and one Cornacchia pass for each modular root found.
    pub fn primitive_representations(
        &self,
        m: &BigUint,
    ) -> Result<Vec<DiagonalBinaryQuadraticRepresentation>, QuadraticFormError> {
        Ok(cornacchia_primitive_solutions(&self.d, m)?
            .into_iter()
            .map(|solution| {
                DiagonalBinaryQuadraticRepresentation::new(
                    solution.x().clone(),
                    solution.y().clone(),
                )
            })
            .collect())
    }

    /// Returns whether this form primitively represents `m`.
    ///
    /// This is deliberately named `primitively_represents` rather than
    /// `represents`: the current implementation answers the primitive
    /// representation problem.
    pub fn primitively_represents(&self, m: &BigUint) -> Result<bool, QuadraticFormError> {
        Ok(!self.primitive_representations(m)?.is_empty())
    }
}
