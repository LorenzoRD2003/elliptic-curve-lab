use crate::elliptic_curves::frobenius::FrobeniusTrace;

/// Exact verification report for the Hasse bound over `F_q`.
///
/// For an elliptic curve over `F_q`, if `t = q + 1 - #E(F_q)` is the trace of
/// the relative Frobenius `π_q`, then Hasse's bound says `|t| <= 2 sqrt(q)`.
///
/// In this implementation, the report stores the underlying [`FrobeniusTrace`]
/// and verifies the equivalent inequality `t^2 <= 4q`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseBoundReport {
    frobenius_trace: FrobeniusTrace,
    trace_square: i128,
    bound_square: i128,
}

impl HasseBoundReport {
    /// Builds a Hasse-bound report from an already computed Frobenius trace.
    ///
    /// Complexity: `Θ(1)`.
    pub fn from_frobenius_trace(frobenius_trace: FrobeniusTrace) -> Self {
        let trace = i128::from(frobenius_trace.trace());
        let trace_square = trace * trace;
        let bound_square = 4 * i128::try_from(frobenius_trace.field_order())
            .expect("stored Frobenius trace should have i128-sized field order");

        Self {
            frobenius_trace,
            trace_square,
            bound_square,
        }
    }

    /// Returns the Frobenius trace package used in the check.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the exact integer value `t^2`.
    pub fn trace_square(&self) -> i128 {
        self.trace_square
    }

    /// Returns the exact integer Hasse bound square `4q`.
    pub fn bound_square(&self) -> i128 {
        self.bound_square
    }

    /// Returns the signed gap `4q - t^2`.
    ///
    /// This value is non-negative exactly when Hasse's inequality holds.
    pub fn slack(&self) -> i128 {
        self.bound_square - self.trace_square
    }

    /// Returns whether the trace satisfies Hasse's bound.
    pub fn holds(&self) -> bool {
        self.trace_square <= self.bound_square
    }
}

impl From<FrobeniusTrace> for HasseBoundReport {
    fn from(value: FrobeniusTrace) -> Self {
        Self::from_frobenius_trace(value)
    }
}
