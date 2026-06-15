use crate::elliptic_curves::{
    CurveError,
    frobenius::{FrobeniusTrace, HasseInterval},
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;

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
    character_sum: i128,
    curve_order: u128,
    trace: i128,
}

impl CharacterSumPointCount {
    /// Builds a validated report from one finite base field and one character sum.
    pub fn new(base_field: FiniteFieldDescriptor, character_sum: i128) -> Result<Self, CurveError> {
        let field_order = field_order_i128(&base_field)?;
        let curve_order_i128 = field_order + 1 + character_sum;
        if curve_order_i128 <= 0 {
            return Err(CurveError::InvalidCurveOrder { order: 0 });
        }

        let curve_order = u128::try_from(curve_order_i128)
            .map_err(|_| CurveError::InvalidCurveOrder { order: 0 })?;
        let trace = -character_sum;

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
    pub fn field_order(&self) -> u128 {
        self.base_field
            .cardinality()
            .expect("stored finite-field descriptor should stay internally consistent")
    }

    /// Returns the character sum `\sum_x χ(f(x))`.
    pub fn character_sum(&self) -> i128 {
        self.character_sum
    }

    /// Returns the resulting curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u128 {
        self.curve_order
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> i128 {
        self.trace
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
        let curve_order = u64::try_from(self.curve_order)
            .map_err(|_| CurveError::InvalidCurveOrder { order: u64::MAX })?;
        FrobeniusTrace::from_order(self.base_field.clone(), curve_order)
    }
}
fn field_order_i128(base_field: &FiniteFieldDescriptor) -> Result<i128, CurveError> {
    base_field
        .cardinality()
        .map_err(|_| CurveError::InvalidFrobeniusBaseField {
            characteristic: base_field.characteristic,
            extension_degree: base_field.extension_degree.get(),
        })
        .and_then(|order| {
            i128::try_from(order).map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: base_field.characteristic,
                extension_degree: base_field.extension_degree.get(),
            })
        })
}
