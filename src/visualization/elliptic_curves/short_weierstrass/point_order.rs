use crate::elliptic_curves::{
    frobenius::group_order::GroupOrderRoute,
    short_weierstrass::point_order::{
        ExhaustivePointOrderReport, HasseIntervalPointOrderReport, PointOrderFromMultipleReport,
        PointOrderReductionStep, PointOrderReport, PointOrderStrategyKind,
    },
};
use crate::visualization::{
    Visualizable,
    elliptic_curves::frobenius::{
        describe_group_order_report, describe_hasse_multiple_search_report,
        format_hasse_multiple_search_report,
    },
};

fn describe_point_order_reduction_step(step: &PointOrderReductionStep) -> String {
    format!(
        "prime {}: exponent in M = {}, removed exponent = {}, remaining multiple = {}",
        step.prime(),
        step.exponent_in_multiple(),
        step.removed_exponent(),
        step.remaining_multiple_after_step()
    )
}

/// Formats an order-from-multiple report compactly.
pub(crate) fn format_point_order_from_multiple_report(
    report: &PointOrderFromMultipleReport,
) -> String {
    format!(
        "ord(P) from M = {} is {}",
        report.supplied_multiple(),
        report.exact_order()
    )
}

/// Describes the prime-peeling order recovery from one known multiple.
pub(crate) fn describe_point_order_from_multiple_report(
    report: &PointOrderFromMultipleReport,
) -> String {
    let mut lines = vec![
        "Point order from multiple".to_string(),
        format!("supplied multiple M: {}", report.supplied_multiple()),
        format!("exact order recovered: {}", report.exact_order()),
        format!("final remaining multiple: {}", report.remaining_multiple()),
        "strategy: divide M by one prime at a time while the smaller multiple still annihilates P"
            .to_string(),
    ];

    for step in report.steps() {
        lines.push(describe_point_order_reduction_step(step));
    }

    lines.join("\n")
}

impl Visualizable for PointOrderFromMultipleReport {
    fn format_compact(&self) -> String {
        format_point_order_from_multiple_report(self)
    }

    fn describe(&self) -> String {
        describe_point_order_from_multiple_report(self)
    }
}

pub(super) fn point_order_strategy_kind_label(strategy: PointOrderStrategyKind) -> &'static str {
    match strategy {
        PointOrderStrategyKind::Exhaustive => "exhaustive",
        PointOrderStrategyKind::FromKnownMultiple => "from known multiple",
        PointOrderStrategyKind::HasseIntervalNaive => "naive Hasse interval",
    }
}

pub(super) fn group_order_route_label_for_order_route(route: GroupOrderRoute) -> &'static str {
    match route {
        GroupOrderRoute::Exhaustive => "exhaustive",
        GroupOrderRoute::QuadraticCharacter => "quadratic character",
        GroupOrderRoute::Schoof => "Schoof",
        GroupOrderRoute::MestreFp => "Mestre",
    }
}

/// Formats an exhaustive point-order report compactly.
fn format_exhaustive_point_order_report(report: &ExhaustivePointOrderReport) -> String {
    format!("ord(P) via exhaustive search = {}", report.exact_order())
}

/// Describes an exhaustive point-order report.
pub(super) fn describe_exhaustive_point_order_report(
    report: &ExhaustivePointOrderReport,
) -> String {
    [
        "Exhaustive point order".to_string(),
        format!("exact order: {}", report.exact_order()),
        "strategy: traverse [n]P in the small ambient finite group until the first identity hit"
            .to_string(),
    ]
    .join("\n")
}

impl Visualizable for ExhaustivePointOrderReport {
    fn format_compact(&self) -> String {
        format_exhaustive_point_order_report(self)
    }

    fn describe(&self) -> String {
        describe_exhaustive_point_order_report(self)
    }
}

/// Formats a Hasse-interval point-order report compactly.
fn format_hasse_interval_point_order_report<P: Visualizable>(
    report: &HasseIntervalPointOrderReport<P>,
) -> String {
    format!(
        "ord(P) via H(q) search = {}",
        report.order_from_multiple().exact_order()
    )
}

/// Describes a Hasse-interval point-order report.
fn describe_hasse_interval_point_order_report<P: Visualizable>(
    report: &HasseIntervalPointOrderReport<P>,
) -> String {
    [
        "Point order via naive Hasse interval".to_string(),
        format!(
            "exact order recovered: {}",
            report.order_from_multiple().exact_order()
        ),
        format!(
            "group-order route: {}",
            group_order_route_label_for_order_route(report.group_order_report().route())
        ),
        describe_group_order_report(report.group_order_report()),
        format!(
            "annihilating-multiple search: {}",
            format_hasse_multiple_search_report(report.multiple_search())
        ),
        describe_hasse_multiple_search_report(report.multiple_search()),
        describe_point_order_from_multiple_report(report.order_from_multiple()),
    ]
    .join("\n")
}

impl<P: Visualizable> Visualizable for HasseIntervalPointOrderReport<P> {
    fn format_compact(&self) -> String {
        format_hasse_interval_point_order_report(self)
    }

    fn describe(&self) -> String {
        describe_hasse_interval_point_order_report(self)
    }
}

/// Formats a unified point-order report compactly.
pub(super) fn format_point_order_report<P: Visualizable>(report: &PointOrderReport<P>) -> String {
    match report {
        PointOrderReport::Exhaustive(report) => format_exhaustive_point_order_report(report),
        PointOrderReport::FromKnownMultiple(report) => {
            format_point_order_from_multiple_report(report)
        }
        PointOrderReport::HasseIntervalNaive(report) => {
            format_hasse_interval_point_order_report(report)
        }
    }
}

/// Describes a unified point-order report.
pub(super) fn describe_point_order_report<P: Visualizable>(report: &PointOrderReport<P>) -> String {
    let mut lines = vec![
        "Point order report".to_string(),
        format!(
            "strategy: {}",
            point_order_strategy_kind_label(report.strategy_kind())
        ),
        format!("exact order: {}", report.exact_order()),
    ];

    match report {
        PointOrderReport::Exhaustive(report) => {
            lines.push(describe_exhaustive_point_order_report(report));
        }
        PointOrderReport::FromKnownMultiple(report) => {
            lines.push(describe_point_order_from_multiple_report(report));
        }
        PointOrderReport::HasseIntervalNaive(report) => {
            lines.push(describe_hasse_interval_point_order_report(report));
        }
    }

    lines.join("\n")
}

impl<P: Visualizable> Visualizable for PointOrderReport<P> {
    fn format_compact(&self) -> String {
        format_point_order_report(self)
    }

    fn describe(&self) -> String {
        describe_point_order_report(self)
    }
}
