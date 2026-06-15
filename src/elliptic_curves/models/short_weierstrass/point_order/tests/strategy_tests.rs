use crate::elliptic_curves::{
    frobenius::group_order::GroupOrderStrategy,
    short_weierstrass::point_order::{PointOrderStrategy, PointOrderStrategyKind},
};

use super::shared::bu;

#[test]
fn strategy_kind_preserves_the_selected_route() {
    assert_eq!(
        PointOrderStrategy::Exhaustive.kind(),
        PointOrderStrategyKind::Exhaustive
    );
    assert_eq!(
        PointOrderStrategy::FromKnownMultiple {
            multiple: bu(6),
            factorization: vec![(bu(2), 1), (bu(3), 1)],
        }
        .kind(),
        PointOrderStrategyKind::FromKnownMultiple
    );
    assert_eq!(
        PointOrderStrategy::HasseIntervalNaive {
            group_order_strategy: GroupOrderStrategy::Auto,
        }
        .kind(),
        PointOrderStrategyKind::HasseIntervalNaive
    );
}
