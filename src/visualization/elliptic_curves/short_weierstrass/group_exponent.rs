use num_bigint::BigUint;

use crate::elliptic_curves::{
    frobenius::group_order::GroupOrderRoute,
    short_weierstrass::group_exponent::{
        ExponentAccumulationReport, ExponentAccumulationStep,
        ExponentLowerBoundGroupOrderVerification, GroupExponentReport, GroupExponentStrategy,
    },
};
use crate::visualization::{
    Visualizable,
    elliptic_curves::{
        frobenius::{describe_group_order_report, format_hasse_interval},
        short_weierstrass::point_order::{
            describe_point_order_report, group_order_route_label_for_order_route,
            point_order_strategy_kind_label,
        },
    },
    shared::yes_no,
};

fn group_exponent_strategy_label(strategy: &GroupExponentStrategy) -> &'static str {
    match strategy {
        GroupExponentStrategy::Exhaustive => "exhaustive",
        GroupExponentStrategy::RandomPoints { .. } => "random points",
    }
}

/// Formats one random-point exponent-accumulation step compactly.
fn format_exponent_accumulation_step<P: Visualizable>(
    step: &ExponentAccumulationStep<P>,
) -> String {
    format!(
        "sample {} gives ord(P) = {}; running lcm = {}",
        step.point().format_compact(),
        step.point_order_report().exact_order(),
        step.accumulated_lcm()
    )
}

/// Describes one random-point exponent-accumulation step.
fn describe_exponent_accumulation_step<P: Visualizable>(
    step: &ExponentAccumulationStep<P>,
) -> String {
    [
        "Exponent accumulation step".to_string(),
        format!("sampled point: {}", step.point().format_compact()),
        format!(
            "point-order route: {}",
            point_order_strategy_kind_label(step.point_order_report().strategy_kind())
        ),
        format!("point order: {}", step.point_order_report().exact_order()),
        format!("running lcm candidate: {}", step.accumulated_lcm()),
        describe_point_order_report(step.point_order_report()),
    ]
    .join("\n")
}

/// Formats the exact exhaustive group-exponent report compactly.
fn format_exhaustive_group_exponent_report(exact_exponent: &BigUint) -> String {
    format!("group exponent = {exact_exponent}")
}

/// Describes the exact exhaustive group-exponent report.
pub(super) fn describe_exhaustive_group_exponent_report(exact_exponent: &BigUint) -> String {
    [
        "Exhaustive group exponent".to_string(),
        format!("exact exponent: {exact_exponent}"),
        "strategy: compute every point order in the tiny ambient group and take their lcm"
            .to_string(),
    ]
    .join("\n")
}

/// Formats the random-point exponent accumulation report compactly.
fn format_exponent_accumulation_report<P: Visualizable>(
    report: &ExponentAccumulationReport<P>,
) -> String {
    format!(
        "group exponent lower bound after {} sample(s) = {}",
        report.samples_taken(),
        report.exponent_lower_bound()
    )
}

/// Describes the random-point exponent accumulation report.
fn describe_exponent_accumulation_report<P: Visualizable>(
    report: &ExponentAccumulationReport<P>,
) -> String {
    let mut lines = vec![
        "Exponent accumulation from random points".to_string(),
        format!(
            "requested samples: {}, completed: {}",
            report.samples_requested(),
            report.samples_taken()
        ),
        format!(
            "completed requested run: {}",
            yes_no(report.completed_requested_samples())
        ),
        format!(
            "point-order route: {}",
            point_order_strategy_kind_label(report.point_order_strategy().kind())
        ),
        format!("exponent lower bound: {}", report.exponent_lower_bound()),
        "interpretation: this lcm is a lower bound for the true exponent and becomes exact only if the sampled orders already capture all prime-power factors"
            .to_string(),
    ];

    for step in report.steps() {
        lines.push(describe_exponent_accumulation_step(step));
    }

    lines.join("\n")
}

impl<P: Visualizable> Visualizable for ExponentAccumulationReport<P> {
    fn format_compact(&self) -> String {
        format_exponent_accumulation_report(self)
    }

    fn describe(&self) -> String {
        describe_exponent_accumulation_report(self)
    }
}

/// Formats an exponent-lower-bound group-order verification compactly.
pub(super) fn format_exponent_lower_bound_group_order_verification(
    report: &ExponentLowerBoundGroupOrderVerification,
) -> String {
    match report.verified_group_order() {
        Some(order) => format!(
            "group order verifies #E(F_q) = {} from lower bound {}",
            order,
            report.exponent_lower_bound()
        ),
        None => format!(
            "group order does not uniquely verify #E(F_q) from lower bound {}",
            report.exponent_lower_bound()
        ),
    }
}

/// Describes an exponent-lower-bound group-order verification.
pub(super) fn describe_exponent_lower_bound_group_order_verification(
    report: &ExponentLowerBoundGroupOrderVerification,
) -> String {
    let mut lines = vec![
        "Exponent lower-bound verification by group order".to_string(),
        format!("exponent lower bound: {}", report.exponent_lower_bound()),
        format!(
            "group-order route: {}",
            group_order_route_label_for_order_route(report.group_order_report().route())
        ),
        format!(
            "Hasse interval: {}",
            format_hasse_interval(&report.group_order_report().hasse_interval())
        ),
        describe_group_order_report(report.group_order_report()),
    ];

    match report.verified_group_order() {
        Some(order) => lines.push(format!(
            "verified group order: {order}\nmeaning: the Hasse interval contains exactly one multiple of the lower bound, so #E(F_q) is forced to equal {order}; this does not by itself certify the exponent"
        )),
        None => lines.push(
            "verified group order: none\nmeaning: the Hasse interval contains zero or at least two multiples of the lower bound, so this check does not force one group order"
                .to_string(),
        ),
    }

    lines.join("\n")
}

impl Visualizable for ExponentLowerBoundGroupOrderVerification {
    fn format_compact(&self) -> String {
        format_exponent_lower_bound_group_order_verification(self)
    }

    fn describe(&self) -> String {
        describe_exponent_lower_bound_group_order_verification(self)
    }
}

/// Formats a unified group-exponent report compactly.
pub(super) fn format_group_exponent_report<P: Visualizable>(
    report: &GroupExponentReport<P>,
) -> String {
    match report {
        GroupExponentReport::Exhaustive(exact_exponent) => {
            format_exhaustive_group_exponent_report(exact_exponent)
        }
        GroupExponentReport::RandomPoints(report) => format_exponent_accumulation_report(report),
    }
}

/// Describes a unified group-exponent report.
pub(super) fn describe_group_exponent_report<P: Visualizable>(
    report: &GroupExponentReport<P>,
) -> String {
    let mut lines = vec![
        "Group exponent report".to_string(),
        format!(
            "strategy: {}",
            group_exponent_strategy_label(&report.strategy())
        ),
        format!("exponent lower bound: {}", report.exponent_lower_bound()),
    ];

    if let Some(exact) = report.exact_exponent() {
        lines.push(format!("exact exponent: {exact}"));
    } else {
        lines.push("exact exponent: not certified by this route".to_string());
    }

    match report {
        GroupExponentReport::Exhaustive(exact_exponent) => {
            lines.push(describe_exhaustive_group_exponent_report(exact_exponent));
        }
        GroupExponentReport::RandomPoints(report) => {
            lines.push(describe_exponent_accumulation_report(report));
        }
    }

    lines.join("\n")
}

impl<P: Visualizable> Visualizable for GroupExponentReport<P> {
    fn format_compact(&self) -> String {
        format_group_exponent_report(self)
    }

    fn describe(&self) -> String {
        describe_group_exponent_report(self)
    }
}
