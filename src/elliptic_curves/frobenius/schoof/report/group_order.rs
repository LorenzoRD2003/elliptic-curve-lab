use num_bigint::BigInt;
use num_traits::ToPrimitive;

use crate::elliptic_curves::{
    CurveError,
    frobenius::schoof::{SchoofTraceCrtOutcome, SchoofTraceCrtReport},
    frobenius::{FrobeniusTrace, HasseInterval},
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use crate::fields::traits::FiniteField;
use crate::numerics::{ceil_div_bigint_by_positive, floor_div_bigint_by_positive};

/// Final outcome of the current Schoof group-order stage.
///
/// This stage starts from one CRT class for the Frobenius trace `t` and then
/// intersects that class with Hasse's trace interval `[-⌊2√q⌋, ⌊2√q⌋]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SchoofGroupOrderOutcome {
    /// The CRT class meets Hasse's trace interval in exactly one integer.
    GroupOrderFound { trace: i128, curve_order: u128 },
    /// The CRT class meets Hasse's trace interval in several integers, so
    /// more odd primes are still needed.
    AmbiguousTraceClass {
        first_trace: i128,
        last_trace: i128,
        candidate_count: u128,
    },
    /// The currently available CRT stage stopped before producing one full
    /// trace class modulo the product of the requested primes.
    BlockedOnOddPrime,
    /// The CRT class did not meet Hasse's trace interval at all.
    ///
    /// Mathematically this should not happen for a correct Schoof run.
    InconsistentWithHasse,
}

/// Route-preserving report for the current end-to-end Schoof group-order stage.
///
/// The stored CRT report remains the source of truth for the prime-by-prime
/// Schoof steps. This wrapper adds only the final Hasse-resolution summary
/// needed to recover `t` and then `#E(F_q) = q + 1 - t`.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofGroupOrderReport<F: FiniteField> {
    base_field: FiniteFieldDescriptor,
    crt_report: SchoofTraceCrtReport<F>,
    outcome: SchoofGroupOrderOutcome,
}

impl<F: FiniteField> SchoofGroupOrderReport<F> {
    pub(crate) fn new(
        base_field: FiniteFieldDescriptor,
        crt_report: SchoofTraceCrtReport<F>,
        outcome: SchoofGroupOrderOutcome,
    ) -> Self {
        Self {
            base_field,
            crt_report,
            outcome,
        }
    }

    /// Returns the finite base-field descriptor.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        self.crt_report.field_order()
    }

    /// Returns the underlying CRT-stage Schoof report.
    pub fn crt_report(&self) -> &SchoofTraceCrtReport<F> {
        &self.crt_report
    }

    /// Returns the discrete Hasse interval `H(q)` for `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn hasse_interval(&self) -> HasseInterval {
        HasseInterval::for_q(self.field_order())
            .expect("stored Schoof field order should define a valid Hasse interval")
    }

    /// Returns the final Hasse-resolution outcome.
    pub fn outcome(&self) -> &SchoofGroupOrderOutcome {
        &self.outcome
    }

    /// Returns the recovered curve order `#E(F_q)` when the current Schoof
    /// data already determines it uniquely.
    pub fn curve_order(&self) -> Option<u128> {
        match self.outcome {
            SchoofGroupOrderOutcome::GroupOrderFound { curve_order, .. } => Some(curve_order),
            _ => None,
        }
    }

    /// Returns the recovered Frobenius trace `t = q + 1 - #E(F_q)` when the
    /// current Schoof data already determines it uniquely.
    pub fn trace(&self) -> Option<i128> {
        match self.outcome {
            SchoofGroupOrderOutcome::GroupOrderFound { trace, .. } => Some(trace),
            _ => None,
        }
    }

    /// Converts the resolved Schoof output into the shared `FrobeniusTrace`
    /// package when the current data already determines one unique trace.
    ///
    /// Complexity: `Θ(1)`.
    pub fn to_frobenius_trace(&self) -> Result<Option<FrobeniusTrace>, CurveError> {
        let Some(curve_order) = self.curve_order() else {
            return Ok(None);
        };
        let curve_order_u64 = u64::try_from(curve_order)
            .map_err(|_| CurveError::InvalidCurveOrder { order: u64::MAX })?;
        FrobeniusTrace::from_order(self.base_field.clone(), curve_order_u64).map(Some)
    }
}

pub(crate) fn finalize_schoof_group_order_report<F: FiniteField>(
    base_field: FiniteFieldDescriptor,
    crt_report: SchoofTraceCrtReport<F>,
) -> Result<SchoofGroupOrderReport<F>, CurveError> {
    let outcome = match crt_report.outcome() {
        SchoofTraceCrtOutcome::BlockedOnOddPrime { .. } => {
            SchoofGroupOrderOutcome::BlockedOnOddPrime
        }
        SchoofTraceCrtOutcome::Combined { solution } => {
            let compatible = hasse_compatible_trace_class(crt_report.field_order(), solution)?;
            match compatible {
                HasseCompatibleTraceClass::Unique { trace } => {
                    let curve_order =
                        crate::elliptic_curves::frobenius::invariants::curve_order_from_field_order_and_trace(
                            crt_report.field_order(),
                            trace,
                        )?;
                    SchoofGroupOrderOutcome::GroupOrderFound { trace, curve_order }
                }
                HasseCompatibleTraceClass::Ambiguous {
                    first_trace,
                    last_trace,
                    candidate_count,
                } => SchoofGroupOrderOutcome::AmbiguousTraceClass {
                    first_trace,
                    last_trace,
                    candidate_count,
                },
                HasseCompatibleTraceClass::Empty => SchoofGroupOrderOutcome::InconsistentWithHasse,
            }
        }
    };

    Ok(SchoofGroupOrderReport::new(base_field, crt_report, outcome))
}

enum HasseCompatibleTraceClass {
    Empty,
    Unique {
        trace: i128,
    },
    Ambiguous {
        first_trace: i128,
        last_trace: i128,
        candidate_count: u128,
    },
}

fn hasse_compatible_trace_class(
    field_order: u128,
    solution: &crate::numerics::chinese_remainder::ChineseRemainderSolution,
) -> Result<HasseCompatibleTraceClass, CurveError> {
    let hasse_interval = HasseInterval::for_q(field_order)?;
    let trace_bound = i128::try_from(hasse_interval.trace_bound())
        .map_err(|_| CurveError::InvalidHasseIntervalFieldOrder { field_order })?;
    let lower = BigInt::from(-trace_bound);
    let upper = BigInt::from(trace_bound);
    let residue = BigInt::from(solution.residue().clone());
    let modulus = BigInt::from(solution.modulus().clone());

    let k_min = ceil_div_bigint_by_positive(&(lower - &residue), &modulus);
    let k_max = floor_div_bigint_by_positive(&(upper - &residue), &modulus);
    if k_min > k_max {
        return Ok(HasseCompatibleTraceClass::Empty);
    }

    let first_trace_bigint = &residue + (&k_min * &modulus);
    let last_trace_bigint = &residue + (&k_max * &modulus);
    let first_trace = bigint_to_i128(&first_trace_bigint, field_order)?;
    let last_trace = bigint_to_i128(&last_trace_bigint, field_order)?;
    let candidate_count = (&k_max - &k_min + BigInt::from(1u8))
        .to_u128()
        .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order })?;

    if candidate_count == 1 {
        Ok(HasseCompatibleTraceClass::Unique { trace: first_trace })
    } else {
        Ok(HasseCompatibleTraceClass::Ambiguous {
            first_trace,
            last_trace,
            candidate_count,
        })
    }
}

fn bigint_to_i128(value: &BigInt, field_order: u128) -> Result<i128, CurveError> {
    value
        .to_i128()
        .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order })
}
