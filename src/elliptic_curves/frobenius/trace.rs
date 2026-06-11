use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::frobenius::HasseInterval;
use crate::fields::FiniteFieldDescriptor;

/// Frobenius trace data recovered from a point count over a finite base field.
///
/// For an elliptic curve over `F_q`, the trace of the relative Frobenius
/// `π_q` is the integer `t = q + 1 - #E(F_q)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusTrace {
    /// Finite base-field metadata for `F_q`.
    base_field: FiniteFieldDescriptor,
    /// The counted order `#E(F_q)`.
    curve_order: u64,
    /// The Frobenius trace `t = q + 1 - #E(F_q)`.
    trace: i64,
}

impl FrobeniusTrace {
    /// Builds a validated Frobenius-trace package from `F_q` and `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`
    pub fn from_order(
        base_field: FiniteFieldDescriptor,
        curve_order: u64,
    ) -> Result<Self, CurveError> {
        let field_order = field_order(&base_field)?;

        if curve_order == 0 {
            return Err(CurveError::InvalidCurveOrder { order: curve_order });
        }

        let trace_i128 = field_order + 1 - i128::from(curve_order);
        let trace = i64::try_from(trace_i128)
            .map_err(|_| CurveError::InvalidCurveOrder { order: curve_order })?;

        Ok(Self {
            base_field,
            curve_order,
            trace,
        })
    }

    /// Reconstructs `#E(F_q)` from `F_q` and the Frobenius trace `t`.
    ///
    /// Complexity: `Θ(1)`
    pub fn curve_order_from_trace(
        base_field: FiniteFieldDescriptor,
        trace: i64,
    ) -> Result<u64, CurveError> {
        let field_order = field_order(&base_field)?;

        let curve_order_i128 = field_order + 1 - i128::from(trace);
        if curve_order_i128 <= 0 {
            return Err(CurveError::InvalidFrobeniusTrace { trace });
        }

        u64::try_from(curve_order_i128).map_err(|_| CurveError::InvalidFrobeniusTrace { trace })
    }

    /// Returns the finite base-field descriptor for `F_q`.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the finite base-field cardinality `q`.
    ///
    /// Complexity: `Θ(1)`
    pub fn field_order(&self) -> u128 {
        self.base_field
            .cardinality()
            .expect("stored finite-field descriptor should stay internally consistent")
    }

    pub fn curve_order(&self) -> u64 {
        self.curve_order
    }

    pub fn trace(&self) -> i64 {
        self.trace
    }

    /// Returns the discrete Hasse interval `H(q)` attached to this base field.
    ///
    /// If this trace package stores data for a curve over `F_q`, the returned
    /// interval is
    ///
    /// `H(q) = [ceil(q + 1 - 2 sqrt(q)), floor(q + 1 + 2 sqrt(q))]`,
    ///
    /// the standard integer search interval that must contain `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn hasse_interval(&self) -> HasseInterval {
        HasseInterval::from_trace(self)
    }
}

fn field_order(base_field: &FiniteFieldDescriptor) -> Result<i128, CurveError> {
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
