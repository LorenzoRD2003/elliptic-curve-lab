use num_bigint::BigUint;

use super::{from_multiple::PointOrderFromMultipleReport, strategies::PointOrderStrategyKind};
use crate::elliptic_curves::frobenius::{
    group_order::GroupOrderReport, hasse::HasseMultipleSearchReport,
};

/// Report for the exhaustive small-group point-order route.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExhaustivePointOrderReport {
    exact_order: BigUint,
}

impl ExhaustivePointOrderReport {
    pub fn new(exact_order: BigUint) -> Self {
        Self { exact_order }
    }

    /// Returns the recovered exact order.
    pub fn exact_order(&self) -> &num_bigint::BigUint {
        &self.exact_order
    }
}

/// Report for the Hasse-interval route to point order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseIntervalPointOrderReport<P> {
    pub(crate) group_order_report: GroupOrderReport,
    pub(crate) multiple_search: HasseMultipleSearchReport<P>,
    pub(crate) order_from_multiple: PointOrderFromMultipleReport,
}

impl<P> HasseIntervalPointOrderReport<P> {
    pub fn group_order_report(&self) -> &GroupOrderReport {
        &self.group_order_report
    }

    pub fn multiple_search(&self) -> &HasseMultipleSearchReport<P> {
        &self.multiple_search
    }

    pub fn order_from_multiple(&self) -> &PointOrderFromMultipleReport {
        &self.order_from_multiple
    }

    pub fn exact_order(&self) -> &num_bigint::BigUint {
        self.order_from_multiple.exact_order()
    }
}

/// Shared point-order report returned by the unified curve-side order API.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointOrderReport<P> {
    Exhaustive(ExhaustivePointOrderReport),
    FromKnownMultiple(PointOrderFromMultipleReport),
    HasseIntervalNaive(Box<HasseIntervalPointOrderReport<P>>),
}

impl<P> PointOrderReport<P> {
    pub fn strategy_kind(&self) -> PointOrderStrategyKind {
        match self {
            Self::Exhaustive(_) => PointOrderStrategyKind::Exhaustive,
            Self::FromKnownMultiple(_) => PointOrderStrategyKind::FromKnownMultiple,
            Self::HasseIntervalNaive(_) => PointOrderStrategyKind::HasseIntervalNaive,
        }
    }

    pub fn exact_order(&self) -> &num_bigint::BigUint {
        match self {
            Self::Exhaustive(report) => report.exact_order(),
            Self::FromKnownMultiple(report) => report.exact_order(),
            Self::HasseIntervalNaive(report) => report.exact_order(),
        }
    }
}
