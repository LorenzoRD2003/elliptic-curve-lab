use crate::elliptic_curves::frobenius::{
    character_sum::CharacterSumPointCount,
    group_order::{
        GroupOrderReport, GroupOrderRoute, MestreGroupOrderReport, MestreSide, MestreStepReport,
    },
    schoof::{
        SchoofGroupOrderOutcome, SchoofGroupOrderReport, SchoofTraceCrtOutcome,
        SchoofTraceCrtReport, SchoofTraceModOddPrimeOutcome,
    },
};
use crate::fields::traits::FiniteField;
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::{
        frobenius::format_hasse_interval,
        short_weierstrass::{
            describe_point_order_from_multiple_report, format_point_order_from_multiple_report,
        },
    },
    polynomials::dense::format_dense_polynomial,
};

/// Formats a character-sum point-count report compactly.
pub(super) fn format_character_sum_point_count(report: &CharacterSumPointCount) -> String {
    format!(
        "#E({}) via χ-sum = {}",
        report.base_field(),
        report.curve_order()
    )
}

/// Describes one quadratic-character point-count report.
pub(super) fn describe_character_sum_point_count(report: &CharacterSumPointCount) -> String {
    [
        "Quadratic-character point count".to_string(),
        format!("base field: {}", report.base_field()),
        format!("field order q: {}", report.field_order()),
        format!("character sum Σ χ(f(x)): {}", report.character_sum()),
        format!("curve order #E(F_q): {}", report.curve_order()),
        format!("trace t = q + 1 - #E(F_q): {}", report.trace()),
        format!(
            "counting formula: #E(F_q) = q + 1 + Σ χ(f(x)) = {} + 1 + {}",
            report.field_order(),
            report.character_sum()
        ),
        "interpretation: this is the Θ(q) character-sum route, distinct from the fully naive affine-pair scan".to_string(),
    ]
    .join("\n")
}

impl Visualizable for CharacterSumPointCount {
    fn format_compact(&self) -> String {
        format_character_sum_point_count(self)
    }

    fn describe(&self) -> String {
        describe_character_sum_point_count(self)
    }
}

fn group_order_route_label(route: GroupOrderRoute) -> &'static str {
    match route {
        GroupOrderRoute::Exhaustive => "exhaustive",
        GroupOrderRoute::QuadraticCharacter => "quadratic character",
        GroupOrderRoute::Schoof => "Schoof",
        GroupOrderRoute::MestreFp => "Mestre",
    }
}

fn mestre_side_label(side: MestreSide) -> &'static str {
    match side {
        MestreSide::Original => "original curve",
        MestreSide::QuadraticTwist => "quadratic twist",
    }
}

fn mestre_iteration_cap_text(report: &MestreGroupOrderReport) -> String {
    match report.config().max_iterations() {
        Some(cap) => format!("iteration cap: {cap}"),
        None => "iteration cap: unbounded".to_string(),
    }
}

/// Formats one alternating Mestre step compactly.
pub(crate) fn format_mestre_step_report(report: &MestreStepReport) -> String {
    format!(
        "{}: M = {}, ord(P) = {}, running λ lower bound = {}",
        mestre_side_label(report.side()),
        report.annihilating_multiple(),
        report.point_order_report().exact_order(),
        report.accumulated_exponent_lower_bound()
    )
}

/// Describes one alternating Mestre step.
pub(crate) fn describe_mestre_step_report(report: &MestreStepReport) -> String {
    [
        "Mestre step".to_string(),
        format!("side: {}", mestre_side_label(report.side())),
        format!(
            "annihilating multiple in H(p): {}",
            report.annihilating_multiple()
        ),
        format!(
            "updated exponent lower bound on this side: {}",
            report.accumulated_exponent_lower_bound()
        ),
        format!(
            "point-order subreport: {}",
            format_point_order_from_multiple_report(report.point_order_report())
        ),
        describe_point_order_from_multiple_report(report.point_order_report()),
    ]
    .join("\n")
}

impl Visualizable for MestreStepReport {
    fn format_compact(&self) -> String {
        format_mestre_step_report(self)
    }

    fn describe(&self) -> String {
        describe_mestre_step_report(self)
    }
}

/// Formats a prime-field Mestre report compactly.
pub(super) fn format_mestre_group_order_report(report: &MestreGroupOrderReport) -> String {
    format!(
        "#E({}) via Mestre = {}",
        report.base_field(),
        report.curve_order()
    )
}

/// Describes a prime-field Mestre report.
pub(super) fn describe_mestre_group_order_report(report: &MestreGroupOrderReport) -> String {
    let mut lines = vec![
        "Mestre group order".to_string(),
        format!("base field: {}", report.base_field()),
        format!("field order p: {}", report.field_order()),
        mestre_iteration_cap_text(report),
        format!("shared Hasse interval: {}", format_hasse_interval(&report.hasse_interval())),
        format!("curve order #E(F_p): {}", report.curve_order()),
        format!("quadratic-twist order #E'(F_p): {}", report.twist_curve_order()),
        format!("trace t = p + 1 - #E(F_p): {}", report.trace()),
        format!("resolved side: {}", report.resolved_side_label()),
        format!(
            "lower bound for λ(E(F_p)): {}",
            report.original_exponent_lower_bound()
        ),
        format!(
            "lower bound for λ(E'(F_p)): {}",
            report.twist_exponent_lower_bound()
        ),
        format!("recorded Mestre steps: {}", report.step_count()),
        "interpretation: Mestre alternates between the curve and one quadratic twist until one side has a unique multiple in H(p)".to_string(),
    ];

    for (index, step) in report.steps().iter().enumerate() {
        lines.push(format!(
            "step {}: {}",
            index + 1,
            format_mestre_step_report(step)
        ));
    }

    lines
        .push("note: the returned group order is always #E(F_p) for the original curve, even if uniqueness was first certified on the twist side".to_string());
    lines.join("\n")
}

impl Visualizable for MestreGroupOrderReport {
    fn format_compact(&self) -> String {
        format_mestre_group_order_report(self)
    }

    fn describe(&self) -> String {
        describe_mestre_group_order_report(self)
    }
}

/// Formats an automatic Schoof group-order summary compactly.
fn format_schoof_group_order_summary(
    report: &crate::elliptic_curves::frobenius::group_order::SchoofGroupOrderSummary,
) -> String {
    format!(
        "#E({}) via Schoof = {}",
        report.resolved().base_field(),
        report.resolved().curve_order()
    )
}

/// Describes an automatic Schoof group-order summary.
fn describe_schoof_group_order_summary(
    report: &crate::elliptic_curves::frobenius::group_order::SchoofGroupOrderSummary,
) -> String {
    [
        "Schoof group order".to_string(),
        format!("base field: {}", report.resolved().base_field()),
        format!("field order q: {}", report.resolved().field_order()),
        format!("curve order #E(F_q): {}", report.resolved().curve_order()),
        format!("trace t = q + 1 - #E(F_q): {}", report.resolved().trace()),
        format!("combined CRT modulus: {}", report.combined_crt_modulus()),
        format!(
            "attempted odd primes: {:?}",
            report.attempted_odd_primes()
        ),
        "interpretation: this summary records the automatic Schoof route that keeps adding odd primes until the CRT modulus forces a unique Hasse-compatible trace".to_string(),
    ]
    .join("\n")
}

/// Formats one detailed Schoof CRT report compactly.
fn format_schoof_trace_crt_report<F: FiniteField>(report: &SchoofTraceCrtReport<F>) -> String {
    match report.outcome() {
        SchoofTraceCrtOutcome::Combined { solution } => {
            format!(
                "Schoof CRT: t ≡ {} (mod {})",
                solution.residue(),
                solution.modulus()
            )
        }
        SchoofTraceCrtOutcome::BlockedOnOddPrime {
            blocked_prime,
            partial_solution,
        } => format!(
            "Schoof CRT blocked at ℓ = {} after t ≡ {} (mod {})",
            blocked_prime,
            partial_solution.residue(),
            partial_solution.modulus()
        ),
    }
}

/// Describes one detailed Schoof CRT report.
fn describe_schoof_trace_crt_report<F>(report: &SchoofTraceCrtReport<F>) -> String
where
    F: FiniteField,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Schoof CRT".to_string(),
        format!("field order q: {}", report.field_order()),
        format!("mod-2 residue: {}", report.mod_2_report().trace_mod_2()),
        format!(
            "resolved congruence count: {}",
            report.resolved_congruences().len()
        ),
    ];

    for odd_prime_report in report.odd_prime_reports() {
        let line = match odd_prime_report.outcome() {
            SchoofTraceModOddPrimeOutcome::TraceFound { trace_mod_ell } => {
                format!(
                    "ℓ = {}: resolved with t ≡ {} (mod ℓ)",
                    odd_prime_report.odd_prime(),
                    trace_mod_ell
                )
            }
            SchoofTraceModOddPrimeOutcome::NonUnitDenominator {
                candidate_trace_mod_ell,
                witness_gcd,
            } => format!(
                "ℓ = {}: skipped after non-unit denominator at candidate {} with gcd witness {}",
                odd_prime_report.odd_prime(),
                candidate_trace_mod_ell,
                format_dense_polynomial(witness_gcd)
            ),
            SchoofTraceModOddPrimeOutcome::ExhaustedCandidates => {
                format!(
                    "ℓ = {}: exhausted candidates without a trace residue",
                    odd_prime_report.odd_prime()
                )
            }
        };
        lines.push(line);
    }

    match report.outcome() {
        SchoofTraceCrtOutcome::Combined { solution } => lines.push(format!(
            "combined CRT class: t ≡ {} (mod {})",
            solution.residue(),
            solution.modulus()
        )),
        SchoofTraceCrtOutcome::BlockedOnOddPrime {
            blocked_prime,
            partial_solution,
        } => lines.push(format!(
            "blocked at ℓ = {} with partial class t ≡ {} (mod {})",
            blocked_prime,
            partial_solution.residue(),
            partial_solution.modulus()
        )),
    }

    lines.join("\n")
}

impl<F: FiniteField> Visualizable for SchoofTraceCrtReport<F>
where
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_schoof_trace_crt_report(self)
    }

    fn describe(&self) -> String {
        describe_schoof_trace_crt_report(self)
    }
}

/// Formats one detailed automatic Schoof group-order report compactly.
fn format_detailed_schoof_group_order_report<F: FiniteField>(
    report: &SchoofGroupOrderReport<F>,
) -> String {
    match report.outcome() {
        SchoofGroupOrderOutcome::GroupOrderFound { curve_order, .. } => {
            format!(
                "#E({}) via automatic Schoof = {}",
                report.base_field(),
                curve_order
            )
        }
        SchoofGroupOrderOutcome::AmbiguousTraceClass {
            candidate_count, ..
        } => format!(
            "automatic Schoof over {} left {} Hasse-compatible trace candidates",
            report.base_field(),
            candidate_count
        ),
        SchoofGroupOrderOutcome::BlockedOnOddPrime => {
            format!(
                "automatic Schoof over {} blocked before resolving the trace",
                report.base_field()
            )
        }
        SchoofGroupOrderOutcome::InconsistentWithHasse => {
            format!(
                "automatic Schoof over {} produced data inconsistent with Hasse",
                report.base_field()
            )
        }
    }
}

/// Describes one detailed automatic Schoof group-order report.
fn describe_detailed_schoof_group_order_report<F: FiniteField>(
    report: &SchoofGroupOrderReport<F>,
) -> String
where
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Automatic Schoof group order".to_string(),
        format!("base field: {}", report.base_field()),
        format!("field order q: {}", report.field_order()),
    ];

    match report.outcome() {
        SchoofGroupOrderOutcome::GroupOrderFound { trace, curve_order } => {
            lines.push(format!("curve order #E(F_q): {}", curve_order));
            lines.push(format!("trace t = q + 1 - #E(F_q): {}", trace));
        }
        SchoofGroupOrderOutcome::AmbiguousTraceClass {
            first_trace,
            last_trace,
            candidate_count,
        } => {
            lines.push(format!(
                "ambiguous trace progression: first = {}, last = {}, count = {}",
                first_trace, last_trace, candidate_count
            ));
        }
        SchoofGroupOrderOutcome::BlockedOnOddPrime => {
            lines.push("outcome: blocked before a full CRT class was assembled".to_string());
        }
        SchoofGroupOrderOutcome::InconsistentWithHasse => {
            lines.push(
                "outcome: the final CRT class is incompatible with Hasse's theorem".to_string(),
            );
        }
    }

    lines.push(describe_schoof_trace_crt_report(report.crt_report()));
    lines.join("\n")
}

impl<F: FiniteField> Visualizable for SchoofGroupOrderReport<F>
where
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_detailed_schoof_group_order_report(self)
    }

    fn describe(&self) -> String {
        describe_detailed_schoof_group_order_report(self)
    }
}

/// Formats a shared point-count report compactly.
pub(super) fn format_group_order_report(report: &GroupOrderReport) -> String {
    match report {
        GroupOrderReport::ExhaustiveTrace(trace) => {
            format!(
                "#E({}) via exhaustive count = {}",
                trace.base_field(),
                trace.curve_order()
            )
        }
        GroupOrderReport::QuadraticCharacter(report) => format_character_sum_point_count(report),
        GroupOrderReport::Schoof(report) => format_schoof_group_order_summary(report),
        GroupOrderReport::MestreFp(report) => format_mestre_group_order_report(report),
    }
}

/// Describes one shared group-order report, including the chosen route.
pub(crate) fn describe_group_order_report(report: &GroupOrderReport) -> String {
    match report {
        GroupOrderReport::ExhaustiveTrace(trace) => [
            "Group order".to_string(),
            format!("strategy: {}", group_order_route_label(report.route())),
            format!("base field: {}", trace.base_field()),
            format!("field order q: {}", trace.field_order()),
            format!("curve order #E(F_q): {}", trace.curve_order()),
            format!("trace t = q + 1 - #E(F_q): {}", trace.trace()),
            "interpretation: this route counts the represented rational points directly by exhaustive enumeration".to_string(),
        ]
        .join("\n"),
        GroupOrderReport::QuadraticCharacter(character_sum) => {
            let mut lines = vec![
                "Group order".to_string(),
                format!("strategy: {}", group_order_route_label(report.route())),
            ];
            lines.extend(
                describe_character_sum_point_count(character_sum)
                    .lines()
                    .skip(1)
                    .map(str::to_string),
            );
            lines.join("\n")
        }
        GroupOrderReport::Schoof(schoof) => {
            let mut lines = vec![
                "Group order".to_string(),
                format!("strategy: {}", group_order_route_label(report.route())),
            ];
            lines.extend(
                describe_schoof_group_order_summary(schoof)
                    .lines()
                    .skip(1)
                    .map(str::to_string),
            );
            lines.join("\n")
        }
        GroupOrderReport::MestreFp(mestre) => {
            let mut lines = vec![
                "Group order".to_string(),
                format!("strategy: {}", group_order_route_label(report.route())),
            ];
            lines.extend(
                describe_mestre_group_order_report(mestre)
                    .lines()
                    .skip(1)
                    .map(str::to_string),
            );
            lines.join("\n")
        }
    }
}

impl Visualizable for GroupOrderReport {
    fn format_compact(&self) -> String {
        format_group_order_report(self)
    }

    fn describe(&self) -> String {
        describe_group_order_report(self)
    }
}
