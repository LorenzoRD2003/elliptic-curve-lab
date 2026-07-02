use crate::elliptic_curves::{
    CurveError,
    frobenius::{FrobeniusTrace, HasseInterval},
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use num_bigint::{BigInt, BigUint};

/// Point-count data recovered from the quadratic-character sum
///
/// `#E(F_q) = q + 1 + Σ_{x ∈ F_q} χ(f(x))`
///
/// for an affine model `y^2 = f(x)` over a finite field of odd characteristic.
///
/// The current implementation is educational and field-enumeration-based:
/// it still loops over all `x ∈ F_q`, but it avoids the fully naive `Θ(q^2)`
/// scan over all affine pairs `(x, y)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CharacterSumPointCount {
    base_field: FiniteFieldDescriptor,
    character_sum: BigInt,
    curve_order: BigUint,
    trace: BigInt,
}

impl CharacterSumPointCount {
    /// Builds a validated report from one finite base field and one character sum.
    pub fn new(
        base_field: FiniteFieldDescriptor,
        character_sum: impl Into<BigInt>,
    ) -> Result<Self, CurveError> {
        let character_sum = character_sum.into();
        let field_order = BigInt::from(base_field.cardinality_biguint());
        let curve_order_bigint = field_order + BigInt::from(1u8) + &character_sum;
        if curve_order_bigint <= BigInt::from(0u8) {
            return Err(CurveError::InvalidCurveOrder {
                order: BigUint::from(0u8),
            });
        }

        let curve_order =
            curve_order_bigint
                .to_biguint()
                .ok_or_else(|| CurveError::InvalidCurveOrder {
                    order: BigUint::from(0u8),
                })?;
        let trace = -&character_sum;

        Ok(Self {
            base_field,
            character_sum,
            curve_order,
            trace,
        })
    }

    /// Returns the finite base-field descriptor.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the field order `q`.
    pub fn field_order(&self) -> BigUint {
        self.base_field.cardinality_biguint()
    }

    /// Returns the character sum `\sum_x χ(f(x))`.
    pub fn character_sum(&self) -> BigInt {
        self.character_sum.clone()
    }

    /// Returns the resulting curve order `#E(F_q)`.
    pub fn curve_order(&self) -> BigUint {
        self.curve_order.clone()
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> BigInt {
        self.trace.clone()
    }

    /// Returns the discrete Hasse interval attached to the same `F_q`.
    pub fn hasse_interval(&self) -> HasseInterval {
        HasseInterval::for_q(self.field_order())
            .expect("stored field order should define a valid Hasse interval")
    }

    /// Converts this count report into the shared Frobenius-trace package.
    ///
    /// This conversion is available only when the counted curve order fits the
    /// current `FrobeniusTrace` representation.
    pub fn to_frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        FrobeniusTrace::from_order(self.base_field.clone(), self.curve_order.clone())
    }
}
