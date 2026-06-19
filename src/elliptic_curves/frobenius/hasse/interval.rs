use std::ops::RangeInclusive;

use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        FrobeniusTrace,
        hasse::{HasseMultipleSearchReport, HasseMultipleSearchStep},
    },
};
use crate::fields::{FieldError, traits::FiniteField};

/// Inputs that can be converted into a finite field order `q` for `H(q)`.
///
/// This small adapter keeps the public constructor ergonomic for both direct
/// integers and checked field-order queries such as `F::order()`.
pub trait HasseFieldOrderInput {
    /// Converts the input into a validated `u128` field order.
    fn into_hasse_field_order(self) -> Result<u128, CurveError>;
}

impl HasseFieldOrderInput for u128 {
    fn into_hasse_field_order(self) -> Result<u128, CurveError> {
        Ok(self)
    }
}

impl HasseFieldOrderInput for Result<u128, FieldError> {
    fn into_hasse_field_order(self) -> Result<u128, CurveError> {
        self.map_err(CurveError::from)
    }
}

/// The discrete Hasse interval of possible values of `#E(F_q)`.
///
/// For an elliptic curve `E/F_q`, Hasse's theorem says that
/// `#E(F_q) = q + 1 - t` with `|t| ≤ 2√q`. Equivalently, since
/// `#E(F_q)` is an integer, `#E(F_q) ∈ H(q)` for the discrete interval
///
/// `H(q) = [ceil(q + 1 - 2 √q), floor(q + 1 + 2 √q)]`.
///
/// The current implementation computes this interval exactly in integer
/// arithmetic via
///
/// `H(q) = [q + 1 - floor(2√q), q + 1 + floor(2√q)]`,
///
/// which avoids floating-point approximations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseInterval {
    q: u128,
    lower: u128,
    upper: u128,
}

impl HasseInterval {
    /// Builds the discrete Hasse interval `H(q)` for a finite field order `q`.
    ///
    /// This constructor rejects values smaller than `2`, since there is no
    /// finite field of order `0` or `1`, and it also rejects arithmetic inputs
    /// that would overflow the exact integer formulas used internally.
    ///
    /// Complexity: `Θ(1)`.
    pub fn for_q(q: impl HasseFieldOrderInput) -> Result<Self, CurveError> {
        let q = q.into_hasse_field_order()?;
        if q < 2 {
            return Err(CurveError::InvalidHasseIntervalFieldOrder { field_order: q });
        }

        let doubled_sqrt_floor = q
            .checked_mul(4)
            .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order: q })?
            .isqrt();
        let center = q
            .checked_add(1)
            .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order: q })?;
        let lower = center
            .checked_sub(doubled_sqrt_floor)
            .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order: q })?;
        let upper = center
            .checked_add(doubled_sqrt_floor)
            .ok_or(CurveError::InvalidHasseIntervalFieldOrder { field_order: q })?;

        Ok(Self { q, lower, upper })
    }

    /// Builds the discrete Hasse interval `H(q)` from one finite field family.
    ///
    /// This is the most direct type-driven entry point when the field itself is
    /// already known at the call site.
    ///
    /// Complexity: `Θ(1)`.
    pub fn for_field<F: FiniteField>() -> Result<Self, CurveError> {
        Self::for_q(F::order())
    }

    /// Rebuilds the interval from an existing Frobenius-trace package.
    ///
    /// Complexity: `Θ(1)`.
    pub fn from_trace(trace: &FrobeniusTrace) -> Self {
        Self::for_q(trace.field_order())
            .expect("stored Frobenius trace should keep the field order valid for H(q)")
    }

    /// Returns the finite field order `q`.
    pub fn q(&self) -> u128 {
        self.q
    }

    /// Returns the lower endpoint of `H(q)`.
    pub fn lower(&self) -> u128 {
        self.lower
    }

    /// Returns the upper endpoint of `H(q)`.
    pub fn upper(&self) -> u128 {
        self.upper
    }

    /// Returns whether `n` lies in the discrete Hasse interval.
    pub fn contains(&self, n: u128) -> bool {
        self.lower <= n && n <= self.upper
    }

    /// Returns `H(q)` as an inclusive integer range.
    pub fn as_range_inclusive(&self) -> RangeInclusive<u128> {
        self.lower..=self.upper
    }

    /// Returns the endpoint difference `upper - lower`.
    ///
    /// This is the discrete span of the stored integer interval, not the
    /// analytic width `4 sqrt(q)`.
    pub fn span(&self) -> u128 {
        self.upper - self.lower
    }

    /// Returns the number of integer candidates inside `H(q)`.
    pub fn candidate_count(&self) -> u128 {
        self.span() + 1
    }

    /// Returns the exact trace bound `⌊2√q⌋` attached to this interval.
    ///
    /// Equivalently, if `#E(F_q) = q + 1 - t` and `#E(F_q) ∈ H(q)`, then the
    /// associated Frobenius trace must satisfy `|t| ≤ trace_bound()`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn trace_bound(&self) -> u128 {
        self.upper - (self.q + 1)
    }

    /// Returns the first multiple of `n` that lies in `H(q)`, if one exists.
    ///
    /// If `n = 0`, this method returns `None`.
    pub fn first_multiple_of(&self, n: u128) -> Option<u128> {
        if n == 0 {
            return None;
        }

        let remainder = self.lower % n;
        let candidate = if remainder == 0 {
            self.lower
        } else {
            self.lower.checked_add(n - remainder)?
        };

        self.contains(candidate).then_some(candidate)
    }

    /// Returns the last multiple of `n` that lies in `H(q)`, if one exists.
    ///
    /// If `n = 0`, this method returns `None`.
    pub fn last_multiple_of(&self, n: u128) -> Option<u128> {
        if n == 0 {
            return None;
        }

        let candidate = self.upper - (self.upper % n);
        self.contains(candidate).then_some(candidate)
    }

    /// Returns how many multiples of `n` lie in `H(q)`.
    ///
    /// If `n = 0`, this method returns `0`.
    pub fn multiple_count_of(&self, n: u128) -> u128 {
        match (self.first_multiple_of(n), self.last_multiple_of(n)) {
            (Some(first), Some(last)) => ((last - first) / n) + 1,
            _ => 0,
        }
    }

    /// Returns the unique multiple of `n` in `H(q)`, if exactly one exists.
    ///
    /// If there are zero or at least two multiples, this method returns `None`.
    pub fn unique_multiple_of(&self, n: u128) -> Option<u128> {
        (self.multiple_count_of(n) == 1)
            .then(|| self.first_multiple_of(n))
            .flatten()
    }

    /// Lists all multiples of `n` that lie in `H(q)`.
    ///
    /// This is mainly an educational helper: algorithms that need only to
    /// distinguish the cases “none”, “unique”, or “several” should prefer
    /// [`Self::multiple_count_of`] or [`Self::unique_multiple_of`].
    ///
    /// If `n = 0`, this method returns the empty vector.
    ///
    /// Complexity: `Θ(k)`, where `k` is the number of returned multiples.
    pub fn multiples_of(&self, n: u128) -> Vec<u128> {
        let Some(first) = self.first_multiple_of(n) else {
            return Vec::new();
        };
        let Some(last) = self.last_multiple_of(n) else {
            return Vec::new();
        };

        let mut multiples = Vec::new();
        let mut current = first;
        loop {
            multiples.push(current);
            if current == last {
                break;
            }
            current = current
                .checked_add(n)
                .expect("iterating between two existing multiples should stay in range");
        }
        multiples
    }

    pub(crate) fn search_report<P>(
        self,
        first_annihilating_multiple: Option<u128>,
        steps: Vec<HasseMultipleSearchStep<P>>,
    ) -> HasseMultipleSearchReport<P> {
        HasseMultipleSearchReport::new(self.q, self, first_annihilating_multiple, steps)
    }
}
