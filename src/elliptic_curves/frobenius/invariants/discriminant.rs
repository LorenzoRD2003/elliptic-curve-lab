use crate::elliptic_curves::{
    endomorphisms::quadratic_orders::QuadraticDiscriminant, frobenius::FrobeniusTrace,
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;

/// Frobenius-side discriminant data derived from one trace package.
///
/// For an elliptic curve over `F_q` with relative Frobenius trace `t`, the
/// Frobenius discriminant is `Δ_π = t^2 - 4q`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusDiscriminant {
    frobenius_trace: FrobeniusTrace,
    discriminant: QuadraticDiscriminant,
}

impl FrobeniusTrace {
    /// Returns the Frobenius discriminant `Δ_π = t^2 - 4q`.
    ///
    /// If this trace package stores the relative Frobenius trace `t` over the
    /// finite base field `F_q`, then the returned report stores the same trace
    /// package together with the derived integral discriminant `t^2 - 4q`.
    ///
    /// Complexity: `Θ(1)` big-integer arithmetic.
    pub fn discriminant(&self) -> FrobeniusDiscriminant {
        FrobeniusDiscriminant::new(self.clone())
    }
}

impl FrobeniusDiscriminant {
    /// Builds the Frobenius discriminant report from one trace package.
    ///
    /// Complexity: `Θ(1)` big-integer arithmetic.
    pub fn new(frobenius_trace: FrobeniusTrace) -> Self {
        let discriminant = QuadraticDiscriminant::from_frobenius_trace_and_field_order(
            frobenius_trace.trace(),
            frobenius_trace.field_order(),
        );

        Self {
            frobenius_trace,
            discriminant,
        }
    }

    /// Returns the originating Frobenius trace package.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the finite base-field descriptor for `F_q`.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        self.frobenius_trace.base_field()
    }

    /// Returns the counted curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u64 {
        self.frobenius_trace.curve_order()
    }

    /// Returns the relative Frobenius trace `t`.
    pub fn trace(&self) -> i64 {
        self.frobenius_trace.trace()
    }

    /// Returns the derived integral quadratic discriminant `Δ_π = t^2 - 4q`.
    pub fn quadratic_discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Returns whether `Δ_π < 0`.
    pub fn is_negative(&self) -> bool {
        self.quadratic_discriminant().is_negative()
    }

    /// Returns whether `Δ_π = 0`.
    pub fn is_zero(&self) -> bool {
        self.quadratic_discriminant().is_zero()
    }

    /// Returns whether `Δ_π > 0`.
    pub fn is_positive(&self) -> bool {
        self.quadratic_discriminant().is_positive()
    }

    /// Returns whether the Frobenius discriminant is fundamental.
    pub fn is_fundamental(&self) -> bool {
        self.quadratic_discriminant().is_fundamental()
    }
}
