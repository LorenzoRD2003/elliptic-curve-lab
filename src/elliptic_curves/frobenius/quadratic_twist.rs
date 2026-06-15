use crate::elliptic_curves::{
    frobenius::FrobeniusTrace, short_weierstrass::isomorphisms::TwistKind,
};

/// Frobenius relation between a curve and a chosen quadratic twist over `F_q`.
///
/// For a non-trivial quadratic twist `E'` of `E` over `F_q`, the classical
/// point-count identity is `#E(F_q) + #E'(F_q) = 2q + 2`.
///
/// Equivalently, if `#E(F_q) = q + 1 - t`, then `#E'(F_q) = q + 1 + t`,
/// so the two Frobenius traces satisfy `t' = -t`.
///
/// This report stores the two Frobenius-trace packages as the primary data and
/// records both the raw invariants and the base-field twist kind attached to
/// the chosen package.
///
/// That distinction matters on the exceptional short-Weierstrass loci
/// `j = 1728` and `j = 0`: a non-square twist factor need not force a
/// genuinely quadratic twist. For example, at `j = 1728` one can have a
/// base-field-trivial twist with a non-square factor `d` whenever `-d` is a
/// square in the base field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticTwistFrobeniusRelation {
    twist_kind: TwistKind,
    original: FrobeniusTrace,
    twist: FrobeniusTrace,
    sum_orders: u128,
    expected_sum: u128,
    holds: bool,
}

impl QuadraticTwistFrobeniusRelation {
    pub(crate) fn new(
        twist_kind: TwistKind,
        original: FrobeniusTrace,
        twist: FrobeniusTrace,
        sum_orders: u128,
        expected_sum: u128,
        holds: bool,
    ) -> Self {
        Self {
            twist_kind,
            original,
            twist,
            sum_orders,
            expected_sum,
            holds,
        }
    }

    /// Returns whether the stored twist package is base-field trivial or
    /// genuinely quadratic.
    pub fn twist_kind(&self) -> TwistKind {
        self.twist_kind
    }

    /// Returns the Frobenius trace package of the original curve `E`.
    pub fn original(&self) -> &FrobeniusTrace {
        &self.original
    }

    /// Returns the Frobenius trace package of the chosen twist `E'`.
    pub fn twist(&self) -> &FrobeniusTrace {
        &self.twist
    }

    /// Returns `#E(F_q) + #E'(F_q)`.
    pub fn sum_orders(&self) -> u128 {
        self.sum_orders
    }

    /// Returns the expected value `2q + 2`.
    pub fn expected_sum(&self) -> u128 {
        self.expected_sum
    }

    /// Returns whether `#E(F_q) + #E'(F_q) = 2q + 2`.
    pub fn holds(&self) -> bool {
        self.holds
    }

    /// Returns whether the two Frobenius traces satisfy `t' = -t`.
    pub fn trace_negation_holds(&self) -> bool {
        self.twist.trace() == -self.original.trace()
    }

    /// Returns whether the two curves have the same point count over `F_q`.
    pub fn same_curve_order_holds(&self) -> bool {
        self.original.curve_order() == self.twist.curve_order()
    }

    /// Returns whether the two Frobenius traces are equal.
    pub fn same_trace_holds(&self) -> bool {
        self.original.trace() == self.twist.trace()
    }

    /// Returns whether the observed Frobenius invariants match the
    /// mathematically expected behavior for the stored twist kind.
    ///
    /// - for `TwistKind::Quadratic`, one expects the classical relation
    ///   `#E(F_q) + #E'(F_q) = 2q + 2`, equivalently `t' = -t`
    /// - for `TwistKind::Trivial`, one expects the two curves to be
    ///   base-field isomorphic, hence to have the same curve order and the
    ///   same trace
    pub fn matches_twist_kind_expectation(&self) -> bool {
        match self.twist_kind {
            TwistKind::Quadratic => self.holds() && self.trace_negation_holds(),
            TwistKind::Trivial => self.same_curve_order_holds() && self.same_trace_holds(),
        }
    }
}
