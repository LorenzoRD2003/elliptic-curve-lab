use num_bigint::BigUint;

use crate::elliptic_curves::frobenius::group_order::SmallFieldGroupOrderStrategy;

/// Public strategy choices for recovering the exact order of one point.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointOrderStrategy {
    Exhaustive,
    FromKnownMultiple {
        multiple: BigUint,
        factorization: Vec<(BigUint, u32)>,
    },
    HasseIntervalNaive {
        group_order_strategy: SmallFieldGroupOrderStrategy,
    },
}

/// Route labels for point-order reports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointOrderStrategyKind {
    Exhaustive,
    FromKnownMultiple,
    HasseIntervalNaive,
}

impl PointOrderStrategy {
    /// Returns the route label without the strategy payload.
    pub fn kind(&self) -> PointOrderStrategyKind {
        match self {
            Self::Exhaustive => PointOrderStrategyKind::Exhaustive,
            Self::FromKnownMultiple { .. } => PointOrderStrategyKind::FromKnownMultiple,
            Self::HasseIntervalNaive { .. } => PointOrderStrategyKind::HasseIntervalNaive,
        }
    }
}
