use std::collections::BTreeMap;

use crate::elliptic_curves::frobenius::character_sum::CharacterSumPointCount;
use crate::elliptic_curves::frobenius::characteristic_equation::{
    FrobeniusCharacteristicEquationCheck, FrobeniusCharacteristicEquationExhaustiveReport,
};
use crate::elliptic_curves::frobenius::extension_counts::{
    FrobeniusExtensionCountReport, FrobeniusExtensionCountSequenceReport,
    FrobeniusExtensionEnumerationComparisonReport,
};
use crate::elliptic_curves::frobenius::group_order::{
    GroupOrderReport, GroupOrderRoute, MestreGroupOrderReport, MestreSide, MestreStepReport,
};
use crate::elliptic_curves::frobenius::hasse::{
    HasseBoundReport, HasseMultipleSearchReport, HasseMultipleSearchStep,
};
use crate::elliptic_curves::frobenius::orbit::FrobeniusOrbit;
use crate::elliptic_curves::frobenius::quadratic_twist::QuadraticTwistFrobeniusRelation;
use crate::elliptic_curves::frobenius::{
    AbsoluteFrobenius, FrobeniusCharacteristicPolynomial, FrobeniusCurveType,
    FrobeniusLocalZetaFunction, FrobeniusTrace, HasseInterval, RelativeFrobenius,
    torsion::{FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport},
};
use crate::elliptic_curves::frobenius::schoof::{
    SchoofGroupOrderOutcome, SchoofGroupOrderReport, SchoofTraceCrtOutcome,
    SchoofTraceCrtReport, SchoofTraceModOddPrimeOutcome,
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use crate::isogenies::frobenius_relation::{
    IsogenyFrobeniusRelation, IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData,
};
use crate::visualization::elliptic_curves::short_weierstrass::{
    describe_point_order_from_multiple_report, format_point_order_from_multiple_report,
};
use crate::visualization::traits::Visualizable;

fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

fn field_symbol(characteristic: u64, extension_degree: u32) -> String {
    if extension_degree == 1 {
        characteristic.to_string()
    } else {
        format!("({characteristic}^{extension_degree})")
    }
}

fn iterated_symbol(base: &str, power: u32) -> String {
    if power == 1 {
        base.to_string()
    } else {
        format!("{base}^{power}")
    }
}

fn orbit_points_text<P: Visualizable>(points: &[P]) -> String {
    if points.is_empty() {
        "[]".to_string()
    } else {
        format!(
            "[{}]",
            points
                .iter()
                .map(Visualizable::format_compact)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

fn descriptor_with_relative_degree(
    base_field: &FiniteFieldDescriptor,
    relative_extension_degree: u32,
) -> String {
    let absolute_extension_degree = base_field
        .extension_degree
        .get()
        .checked_mul(relative_extension_degree)
        .expect("small visualization examples should keep extension degrees in range");
    if absolute_extension_degree == 1 {
        format!("F_{}", base_field.characteristic)
    } else {
        format!(
            "F_({}^{absolute_extension_degree})",
            base_field.characteristic
        )
    }
}

#[cfg_attr(not(test), allow(dead_code))]
fn point_power_histogram<P: PartialEq>(
    report: &FrobeniusOnExactTorsionReport<P>,
) -> Option<String> {
    let mut histogram = BTreeMap::<u32, usize>::new();
    for point in report.points() {
        let minimal_power = point.minimal_absolute_frobenius_fixing_power()?;
        *histogram.entry(minimal_power).or_default() += 1;
    }

    Some(
        histogram
            .into_iter()
            .map(|(power, count)| format!("d = {power}: {count}"))
            .collect::<Vec<_>>()
            .join(", "),
    )
}

/// Formats the absolute Frobenius metadata compactly.
pub fn format_absolute_frobenius(frobenius: &AbsoluteFrobenius) -> String {
    iterated_symbol(
        &format!("π_{}", frobenius.characteristic()),
        frobenius.power(),
    )
}

/// Describes the absolute Frobenius metadata.
pub fn describe_absolute_frobenius(frobenius: &AbsoluteFrobenius) -> String {
    [
        "Absolute Frobenius".to_string(),
        format!("symbol: {}", format_absolute_frobenius(frobenius)),
        format!("characteristic p: {}", frobenius.characteristic()),
        format!("iterate: {}", frobenius.power()),
        format!("identity iterate: {}", yes_no(frobenius.is_identity())),
        "interpretation: this is the coordinate-wise p-power map on the represented finite-field backend".to_string(),
    ]
    .join("\n")
}

impl Visualizable for AbsoluteFrobenius {
    fn format_compact(&self) -> String {
        format_absolute_frobenius(self)
    }

    fn describe(&self) -> String {
        describe_absolute_frobenius(self)
    }
}

/// Formats the relative Frobenius metadata compactly.
pub fn format_relative_frobenius(frobenius: &RelativeFrobenius) -> String {
    let symbol = field_symbol(
        frobenius.base_field().characteristic,
        frobenius.base_field().extension_degree.get(),
    );
    iterated_symbol(&format!("π_{symbol}"), frobenius.power())
}

/// Describes the relative Frobenius metadata.
pub fn describe_relative_frobenius(frobenius: &RelativeFrobenius) -> String {
    [
        "Relative Frobenius".to_string(),
        format!("symbol: {}", format_relative_frobenius(frobenius)),
        format!("base field: {}", frobenius.base_field()),
        format!(
            "field order q: {}",
            frobenius
                .base_field()
                .cardinality()
                .expect("stored field descriptors should stay consistent")
        ),
        format!("iterate: {}", frobenius.power()),
        format!("identity iterate: {}", yes_no(frobenius.is_identity())),
        "interpretation: this is the q-power Frobenius attached to the chosen represented base field".to_string(),
    ]
    .join("\n")
}

impl Visualizable for RelativeFrobenius {
    fn format_compact(&self) -> String {
        format_relative_frobenius(self)
    }

    fn describe(&self) -> String {
        describe_relative_frobenius(self)
    }
}

/// Formats Frobenius trace data compactly.
pub fn format_frobenius_trace(trace: &FrobeniusTrace) -> String {
    format!("t = {} over {}", trace.trace(), trace.base_field())
}

/// Describes a Frobenius trace package recovered from `#E(F_q)`.
pub fn describe_frobenius_trace(trace: &FrobeniusTrace) -> String {
    [
        "Frobenius trace".to_string(),
        format!("base field: {}", trace.base_field()),
        format!("field order q: {}", trace.field_order()),
        format!("curve order #E(F_q): {}", trace.curve_order()),
        format!("trace t: {}", trace.trace()),
        format!(
            "defining formula: t = q + 1 - #E(F_q) = {} + 1 - {}",
            trace.field_order(),
            trace.curve_order()
        ),
    ]
    .join("\n")
}

impl Visualizable for FrobeniusTrace {
    fn format_compact(&self) -> String {
        format_frobenius_trace(self)
    }

    fn describe(&self) -> String {
        describe_frobenius_trace(self)
    }
}

/// Formats a Hasse interval compactly.
pub fn format_hasse_interval(interval: &HasseInterval) -> String {
    format!(
        "H({}) = [{} , {}]",
        interval.q(),
        interval.lower(),
        interval.upper()
    )
}

/// Describes the discrete Hasse interval of possible values of `#E(F_q)`.
pub fn describe_hasse_interval(interval: &HasseInterval) -> String {
    let doubled_sqrt_floor = interval.upper() - (interval.q() + 1);
    [
        "Hasse interval".to_string(),
        format!("field order q: {}", interval.q()),
        format!("interval H(q): [{} , {}]", interval.lower(), interval.upper()),
        format!("endpoint gap upper - lower: {}", interval.span()),
        format!("integer candidate count: {}", interval.candidate_count()),
        format!("floor(sqrt(4q)): {}", doubled_sqrt_floor),
        format!(
            "integer formula: H(q) = [q + 1 - floor(sqrt(4q)), q + 1 + floor(sqrt(4q))] = [{} , {}]",
            interval.lower(),
            interval.upper()
        ),
        "interpretation: every possible value of #E(F_q) lies in this discrete interval by Hasse's theorem".to_string(),
    ]
    .join("\n")
}

impl Visualizable for HasseInterval {
    fn format_compact(&self) -> String {
        format_hasse_interval(self)
    }

    fn describe(&self) -> String {
        describe_hasse_interval(self)
    }
}

/// Formats one tested candidate in a Hasse-interval multiple search.
pub fn format_hasse_multiple_search_step<P: Visualizable>(
    step: &HasseMultipleSearchStep<P>,
) -> String {
    format!(
        "M = {} gives [M]P = {}",
        step.candidate_multiple(),
        step.image().format_compact()
    )
}

/// Describes one tested candidate in a Hasse-interval multiple search.
pub fn describe_hasse_multiple_search_step<P: Visualizable>(
    step: &HasseMultipleSearchStep<P>,
) -> String {
    [
        "Hasse-interval multiple search step".to_string(),
        format!("candidate multiple M: {}", step.candidate_multiple()),
        format!("image [M]P: {}", step.image().format_compact()),
    ]
    .join("\n")
}

/// Formats a Hasse-interval annihilating-multiple search compactly.
pub fn format_hasse_multiple_search_report<P: Visualizable>(
    report: &HasseMultipleSearchReport<P>,
) -> String {
    match report.first_annihilating_multiple() {
        Some(multiple) => format!("first H(q)-multiple annihilating P: {multiple}"),
        None => format!(
            "no annihilating multiple found inside {}",
            report.interval().format_compact()
        ),
    }
}

/// Describes a Hasse-interval annihilating-multiple search.
pub fn describe_hasse_multiple_search_report<P: Visualizable>(
    report: &HasseMultipleSearchReport<P>,
) -> String {
    let mut lines = vec![
        "Hasse-interval annihilating-multiple search".to_string(),
        format!("field order q: {}", report.q()),
        format!("searched interval: {}", report.interval().format_compact()),
        format!("tested candidates: {}", report.tested_candidates()),
    ];

    match report.first_annihilating_multiple() {
        Some(multiple) => lines.push(format!("first annihilating multiple: {multiple}")),
        None => lines.push("first annihilating multiple: none found".to_string()),
    }

    lines.push("tested steps:".to_string());
    lines.extend(
        report
            .steps()
            .iter()
            .map(|step| format!("  {}", format_hasse_multiple_search_step(step))),
    );

    lines.join("\n")
}

impl<P: Visualizable> Visualizable for HasseMultipleSearchStep<P> {
    fn format_compact(&self) -> String {
        format_hasse_multiple_search_step(self)
    }

    fn describe(&self) -> String {
        describe_hasse_multiple_search_step(self)
    }
}

impl<P: Visualizable> Visualizable for HasseMultipleSearchReport<P> {
    fn format_compact(&self) -> String {
        format_hasse_multiple_search_report(self)
    }

    fn describe(&self) -> String {
        describe_hasse_multiple_search_report(self)
    }
}

/// Formats a character-sum point-count report compactly.
pub fn format_character_sum_point_count(report: &CharacterSumPointCount) -> String {
    format!(
        "#E({}) via χ-sum = {}",
        report.base_field(),
        report.curve_order()
    )
}

/// Describes one quadratic-character point-count report.
pub fn describe_character_sum_point_count(report: &CharacterSumPointCount) -> String {
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
pub fn format_mestre_group_order_report(report: &MestreGroupOrderReport) -> String {
    format!(
        "#E({}) via Mestre = {}",
        report.base_field(),
        report.curve_order()
    )
}

/// Describes a prime-field Mestre report.
pub fn describe_mestre_group_order_report(report: &MestreGroupOrderReport) -> String {
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
pub fn format_schoof_group_order_summary(
    report: &crate::elliptic_curves::frobenius::group_order::SchoofGroupOrderSummary,
) -> String {
    format!(
        "#E({}) via Schoof = {}",
        report.resolved().base_field(),
        report.resolved().curve_order()
    )
}

/// Describes an automatic Schoof group-order summary.
pub fn describe_schoof_group_order_summary(
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
pub fn format_schoof_trace_crt_report<F: crate::fields::traits::FiniteField>(
    report: &SchoofTraceCrtReport<F>,
) -> String {
    match report.outcome() {
        SchoofTraceCrtOutcome::Combined { solution } => {
            format!("Schoof CRT: t ≡ {} (mod {})", solution.residue(), solution.modulus())
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
pub fn describe_schoof_trace_crt_report<F: crate::fields::traits::FiniteField>(
    report: &SchoofTraceCrtReport<F>,
) -> String {
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
                format!("ℓ = {}: resolved with t ≡ {} (mod ℓ)", odd_prime_report.odd_prime(), trace_mod_ell)
            }
            SchoofTraceModOddPrimeOutcome::NonUnitDenominator {
                candidate_trace_mod_ell,
                witness_gcd,
            } => format!(
                "ℓ = {}: skipped after non-unit denominator at candidate {} with gcd witness degree {:?}",
                odd_prime_report.odd_prime(),
                candidate_trace_mod_ell,
                witness_gcd.degree()
            ),
            SchoofTraceModOddPrimeOutcome::ExhaustedCandidates => {
                format!("ℓ = {}: exhausted candidates without a trace residue", odd_prime_report.odd_prime())
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

impl<F: crate::fields::traits::FiniteField> Visualizable for SchoofTraceCrtReport<F> {
    fn format_compact(&self) -> String {
        format_schoof_trace_crt_report(self)
    }

    fn describe(&self) -> String {
        describe_schoof_trace_crt_report(self)
    }
}

/// Formats one detailed automatic Schoof group-order report compactly.
pub fn format_detailed_schoof_group_order_report<F: crate::fields::traits::FiniteField>(
    report: &SchoofGroupOrderReport<F>,
) -> String {
    match report.outcome() {
        SchoofGroupOrderOutcome::GroupOrderFound { curve_order, .. } => {
            format!("#E({}) via automatic Schoof = {}", report.base_field(), curve_order)
        }
        SchoofGroupOrderOutcome::AmbiguousTraceClass { candidate_count, .. } => format!(
            "automatic Schoof over {} left {} Hasse-compatible trace candidates",
            report.base_field(),
            candidate_count
        ),
        SchoofGroupOrderOutcome::BlockedOnOddPrime => {
            format!("automatic Schoof over {} blocked before resolving the trace", report.base_field())
        }
        SchoofGroupOrderOutcome::InconsistentWithHasse => {
            format!("automatic Schoof over {} produced data inconsistent with Hasse", report.base_field())
        }
    }
}

/// Describes one detailed automatic Schoof group-order report.
pub fn describe_detailed_schoof_group_order_report<F: crate::fields::traits::FiniteField>(
    report: &SchoofGroupOrderReport<F>,
) -> String {
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
            lines.push("outcome: the final CRT class is incompatible with Hasse's theorem".to_string());
        }
    }

    lines.push(describe_schoof_trace_crt_report(report.crt_report()));
    lines.join("\n")
}

impl<F: crate::fields::traits::FiniteField> Visualizable for SchoofGroupOrderReport<F> {
    fn format_compact(&self) -> String {
        format_detailed_schoof_group_order_report(self)
    }

    fn describe(&self) -> String {
        describe_detailed_schoof_group_order_report(self)
    }
}

/// Formats a shared point-count report compactly.
pub fn format_group_order_report(report: &GroupOrderReport) -> String {
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
pub fn describe_group_order_report(report: &GroupOrderReport) -> String {
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

/// Describes the Frobenius characteristic polynomial.
pub fn describe_frobenius_characteristic_polynomial(
    polynomial: &FrobeniusCharacteristicPolynomial,
) -> String {
    [
        "Frobenius characteristic polynomial".to_string(),
        format!("base field: {}", polynomial.base_field()),
        format!("field order q: {}", polynomial.field_order()),
        format!("trace t: {}", polynomial.trace()),
        format!("polynomial: χ_π(T) = {}", polynomial.pretty()),
        format!(
            "evaluation at T = 1: χ_π(1) = {}",
            polynomial.evaluate_at_integer(1)
        ),
        format!("discriminant t^2 - 4q: {}", polynomial.discriminant()),
    ]
    .join("\n")
}

impl Visualizable for FrobeniusCharacteristicPolynomial {
    fn format_compact(&self) -> String {
        self.pretty()
    }

    fn describe(&self) -> String {
        describe_frobenius_characteristic_polynomial(self)
    }
}

/// Describes the local zeta function attached to a Frobenius characteristic polynomial.
pub fn describe_frobenius_local_zeta_function(zeta: &FrobeniusLocalZetaFunction) -> String {
    [
        "Local zeta function".to_string(),
        format!("base field: {}", zeta.base_field()),
        format!("field order q: {}", zeta.field_order()),
        format!("trace t: {}", zeta.trace()),
        format!(
            "source characteristic polynomial: {}",
            zeta.characteristic_polynomial()
        ),
        format!("numerator: {}", zeta.numerator_pretty()),
        format!("denominator: {}", zeta.denominator_pretty()),
        format!("zeta function: {}", zeta.pretty()),
    ]
    .join("\n")
}

impl Visualizable for FrobeniusLocalZetaFunction {
    fn format_compact(&self) -> String {
        self.pretty()
    }

    fn describe(&self) -> String {
        describe_frobenius_local_zeta_function(self)
    }
}

fn curve_type_text(curve_type: FrobeniusCurveType) -> &'static str {
    match curve_type {
        FrobeniusCurveType::Ordinary => "ordinary",
        FrobeniusCurveType::Supersingular => "supersingular",
    }
}

impl Visualizable for FrobeniusCurveType {
    fn format_compact(&self) -> String {
        curve_type_text(*self).to_string()
    }

    fn describe(&self) -> String {
        match self {
            Self::Ordinary => {
                "Ordinary curve: the base-field characteristic does not divide the Frobenius trace."
                    .to_string()
            }
            Self::Supersingular => {
                "Supersingular curve: the base-field characteristic divides the Frobenius trace."
                    .to_string()
            }
        }
    }
}

/// Describes an exact Hasse-bound report.
pub fn describe_hasse_bound_report(report: &HasseBoundReport) -> String {
    [
        "Hasse bound".to_string(),
        format!("base field: {}", report.frobenius_trace().base_field()),
        format!("field order q: {}", report.frobenius_trace().field_order()),
        format!("trace t: {}", report.frobenius_trace().trace()),
        format!("trace square t^2: {}", report.trace_square()),
        format!("bound square 4q: {}", report.bound_square()),
        format!("slack 4q - t^2: {}", report.slack()),
        format!("holds: {}", yes_no(report.holds())),
    ]
    .join("\n")
}

impl Visualizable for HasseBoundReport {
    fn format_compact(&self) -> String {
        format!(
            "Hasse over {}: {}",
            self.frobenius_trace().base_field(),
            yes_no(self.holds())
        )
    }

    fn describe(&self) -> String {
        describe_hasse_bound_report(self)
    }
}

/// Formats one Frobenius orbit compactly.
pub fn format_frobenius_orbit<P: Visualizable>(orbit: &FrobeniusOrbit<P>) -> String {
    format!(
        "period {} orbit {}",
        orbit.period(),
        orbit_points_text(orbit.points())
    )
}

/// Describes one Frobenius orbit.
pub fn describe_frobenius_orbit<P: Visualizable>(orbit: &FrobeniusOrbit<P>) -> String {
    [
        "Frobenius orbit".to_string(),
        format!("start: {}", orbit.start().format_compact()),
        format!("period: {}", orbit.period()),
        format!("points: {}", orbit_points_text(orbit.points())),
    ]
    .join("\n")
}

impl<P: Visualizable> Visualizable for FrobeniusOrbit<P> {
    fn format_compact(&self) -> String {
        format_frobenius_orbit(self)
    }

    fn describe(&self) -> String {
        describe_frobenius_orbit(self)
    }
}

/// Describes one pointwise Frobenius-on-torsion datum.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn describe_frobenius_on_exact_torsion_point<P: Visualizable + PartialEq>(
    point: &FrobeniusOnExactTorsionPoint<P>,
) -> String {
    let mut lines = vec![
        "Frobenius on exact torsion point".to_string(),
        format!("point: {}", point.point().format_compact()),
        format!("image: {}", point.frobenius_image().format_compact()),
        format!(
            "fixed by chosen Frobenius: {}",
            yes_no(point.fixed_by_frobenius())
        ),
    ];

    match point.minimal_absolute_frobenius_fixing_power() {
        Some(power) => lines.push(format!(
            "minimal absolute-Frobenius fixing power: {}",
            power
        )),
        None => lines.push(
            "minimal absolute-Frobenius fixing power: not recorded on this report".to_string(),
        ),
    }

    lines.join("\n")
}

impl<P: Visualizable + PartialEq> Visualizable for FrobeniusOnExactTorsionPoint<P> {
    fn format_compact(&self) -> String {
        format!(
            "{} -> {} ({})",
            self.point().format_compact(),
            self.frobenius_image().format_compact(),
            if self.fixed_by_frobenius() {
                "fixed"
            } else {
                "moved"
            }
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_on_exact_torsion_point(self)
    }
}

/// Describes a Frobenius action on exact rational torsion.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn describe_frobenius_on_exact_torsion_report<P>(
    report: &FrobeniusOnExactTorsionReport<P>,
) -> String
where
    P: Visualizable + Clone + Eq + std::hash::Hash + PartialEq,
{
    let mut lines = vec![
        "Frobenius on exact torsion".to_string(),
        format!("exact order n: {}", report.exact_order()),
        format!("listed torsion points: {}", report.points().len()),
        format!("all fixed: {}", yes_no(report.all_fixed())),
        format!("fixed count: {}", report.fixed_count()),
        format!("moved count: {}", report.moved_count()),
        format!("orbit count: {}", report.orbit_count()),
        format!("orbit periods: {:?}", report.orbit_periods()),
    ];

    match point_power_histogram(report) {
        Some(histogram) => lines.push(format!(
            "minimal absolute-Frobenius fixing powers: {histogram}"
        )),
        None => lines.push(
            "minimal absolute-Frobenius fixing powers: not recorded on this relative-Frobenius report"
                .to_string(),
        ),
    }

    lines.push(format!(
        "pointwise data: {}",
        report
            .points()
            .iter()
            .map(Visualizable::format_compact)
            .collect::<Vec<_>>()
            .join(", ")
    ));

    lines.join("\n")
}

impl<P> Visualizable for FrobeniusOnExactTorsionReport<P>
where
    P: Visualizable + Clone + Eq + std::hash::Hash + PartialEq,
{
    fn format_compact(&self) -> String {
        format!(
            "exact-{} torsion under Frobenius: {} points, {} fixed",
            self.exact_order(),
            self.points().len(),
            self.fixed_count()
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_on_exact_torsion_report(self)
    }
}

/// Describes one point-count report over `F_{q^n}` derived from Frobenius data.
pub fn describe_frobenius_extension_count_report(report: &FrobeniusExtensionCountReport) -> String {
    let extension_field = descriptor_with_relative_degree(
        report.frobenius_trace().base_field(),
        report.extension_degree().get(),
    );
    [
        "Frobenius extension count".to_string(),
        format!("base field: {}", report.frobenius_trace().base_field()),
        format!("extension field: {extension_field}"),
        format!("relative extension degree n: {}", report.extension_degree()),
        format!(
            "extension field order q^n: {}",
            report.extension_field_order()
        ),
        format!("power sum s_n = α^n + β^n: {}", report.power_sum()),
        format!(
            "curve order #E({extension_field}): {}",
            report.curve_order()
        ),
    ]
    .join("\n")
}

impl Visualizable for FrobeniusExtensionCountReport {
    fn format_compact(&self) -> String {
        let extension_field = descriptor_with_relative_degree(
            self.frobenius_trace().base_field(),
            self.extension_degree().get(),
        );
        format!("#E({extension_field}) = {}", self.curve_order())
    }

    fn describe(&self) -> String {
        describe_frobenius_extension_count_report(self)
    }
}

/// Describes a prefix of extension counts derived from one Frobenius trace.
pub fn describe_frobenius_extension_count_sequence_report(
    report: &FrobeniusExtensionCountSequenceReport,
) -> String {
    let mut lines = vec![
        "Frobenius extension count sequence".to_string(),
        format!("base field: {}", report.frobenius_trace().base_field()),
        format!("stored degrees: {}", report.reports().len()),
    ];

    for count in report.reports() {
        let extension_field = descriptor_with_relative_degree(
            report.frobenius_trace().base_field(),
            count.extension_degree().get(),
        );
        lines.push(format!(
            "degree {}: #E({extension_field}) = {}",
            count.extension_degree(),
            count.curve_order()
        ));
    }

    lines.join("\n")
}

impl Visualizable for FrobeniusExtensionCountSequenceReport {
    fn format_compact(&self) -> String {
        format!(
            "extension counts through n = {} over {}",
            self.reports().len(),
            self.frobenius_trace().base_field()
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_extension_count_sequence_report(self)
    }
}

/// Describes the comparison between the Frobenius-derived and exhaustive extension counts.
pub fn describe_frobenius_extension_enumeration_comparison_report(
    report: &FrobeniusExtensionEnumerationComparisonReport,
) -> String {
    [
        "Frobenius versus exhaustive extension count".to_string(),
        format!("trace base field: {}", report.trace_base_field()),
        format!("represented curve base field: {}", report.curve_base_field()),
        format!(
            "relative extension degree: {}",
            report.relative_extension_degree()
        ),
        format!(
            "Frobenius-derived count: {}",
            report.frobenius_count().curve_order()
        ),
        format!("exhaustive count: {}", report.exhaustive_curve_order()),
        format!("agreement: {}", yes_no(report.agrees())),
        "interpretation: this report keeps the fast Frobenius recurrence and the slow exhaustive enumeration visible side by side".to_string(),
    ]
    .join("\n")
}

impl Visualizable for FrobeniusExtensionEnumerationComparisonReport {
    fn format_compact(&self) -> String {
        format!(
            "Frobenius vs exhaustive over {}: {}",
            self.curve_base_field(),
            yes_no(self.agrees())
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_extension_enumeration_comparison_report(self)
    }
}

/// Describes one pointwise check of the Frobenius characteristic equation.
pub fn describe_frobenius_characteristic_equation_check<P: Visualizable>(
    check: &FrobeniusCharacteristicEquationCheck<P>,
) -> String {
    [
        "Frobenius characteristic equation at one point".to_string(),
        format!("point P: {}", check.point().format_compact()),
        format!("π_q(P): {}", check.pi_q().format_compact()),
        format!("π_q^2(P): {}", check.pi_q_squared().format_compact()),
        format!("[t]π_q(P): {}", check.trace_term().format_compact()),
        format!("[q]P: {}", check.q_times_point().format_compact()),
        format!("left-hand side: {}", check.lhs().format_compact()),
        format!("holds: {}", yes_no(check.holds())),
    ]
    .join("\n")
}

impl<P: Visualizable> Visualizable for FrobeniusCharacteristicEquationCheck<P> {
    fn format_compact(&self) -> String {
        format!(
            "characteristic equation at {}: {}",
            self.point().format_compact(),
            yes_no(self.holds())
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_characteristic_equation_check(self)
    }
}

/// Describes an exhaustive characteristic-equation verification report.
pub fn describe_frobenius_characteristic_equation_exhaustive_report<P: Visualizable>(
    report: &FrobeniusCharacteristicEquationExhaustiveReport<P>,
) -> String {
    let failed_points = if report.failed_checks().is_empty() {
        "none".to_string()
    } else {
        report
            .failed_checks()
            .iter()
            .map(|check| check.point().format_compact())
            .collect::<Vec<_>>()
            .join(", ")
    };

    [
        "Exhaustive Frobenius characteristic equation check".to_string(),
        format!("base field: {}", report.frobenius_trace().base_field()),
        format!("trace t: {}", report.frobenius_trace().trace()),
        format!("checked points: {}", report.checked_points()),
        format!("all hold: {}", yes_no(report.all_hold())),
        format!("failed checks: {}", report.failed_checks().len()),
        format!("failed points: {failed_points}"),
    ]
    .join("\n")
}

impl<P: Visualizable> Visualizable for FrobeniusCharacteristicEquationExhaustiveReport<P> {
    fn format_compact(&self) -> String {
        format!(
            "exhaustive characteristic equation over {}: {}",
            self.frobenius_trace().base_field(),
            yes_no(self.all_hold())
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_characteristic_equation_exhaustive_report(self)
    }
}

/// Describes the Frobenius relation between a curve and a chosen quadratic twist.
pub fn describe_quadratic_twist_frobenius_relation(
    relation: &QuadraticTwistFrobeniusRelation,
) -> String {
    [
        "Quadratic-twist Frobenius relation".to_string(),
        format!("base field: {}", relation.original().base_field()),
        format!(
            "twist kind over the base field: {:?}",
            relation.twist_kind()
        ),
        format!("original order: {}", relation.original().curve_order()),
        format!("twist order: {}", relation.twist().curve_order()),
        format!("original trace: {}", relation.original().trace()),
        format!("twist trace: {}", relation.twist().trace()),
        format!("sum of orders: {}", relation.sum_orders()),
        format!("expected sum 2q + 2: {}", relation.expected_sum()),
        format!("order relation holds: {}", yes_no(relation.holds())),
        format!(
            "same curve order holds: {}",
            yes_no(relation.same_curve_order_holds())
        ),
        format!("same trace holds: {}", yes_no(relation.same_trace_holds())),
        format!(
            "trace negation t' = -t holds: {}",
            yes_no(relation.trace_negation_holds())
        ),
        format!(
            "matches the expectation for this twist kind: {}",
            yes_no(relation.matches_twist_kind_expectation())
        ),
    ]
    .join("\n")
}

impl Visualizable for QuadraticTwistFrobeniusRelation {
    fn format_compact(&self) -> String {
        format!(
            "quadratic-twist Frobenius relation over {} ({:?}): {}",
            self.original().base_field(),
            self.twist_kind(),
            yes_no(self.matches_twist_kind_expectation())
        )
    }

    fn describe(&self) -> String {
        describe_quadratic_twist_frobenius_relation(self)
    }
}

/// Describes the Frobenius relation attached to an explicit isogeny.
pub fn describe_isogeny_frobenius_relation(relation: &IsogenyFrobeniusRelation) -> String {
    [
        "Isogeny Frobenius relation".to_string(),
        format!("base field: {}", relation.domain().base_field()),
        format!("isogeny degree: {}", relation.degree()),
        format!("domain order: {}", relation.domain().curve_order()),
        format!("codomain order: {}", relation.codomain().curve_order()),
        format!("domain trace: {}", relation.domain().trace()),
        format!("codomain trace: {}", relation.codomain().trace()),
        format!("same curve order: {}", yes_no(relation.same_curve_order())),
        format!("same trace: {}", yes_no(relation.same_trace())),
        format!("holds: {}", yes_no(relation.holds())),
    ]
    .join("\n")
}

impl Visualizable for IsogenyFrobeniusRelation {
    fn format_compact(&self) -> String {
        format!(
            "isogeny degree {} over {}: {}",
            self.degree(),
            self.domain().base_field(),
            yes_no(self.holds())
        )
    }

    fn describe(&self) -> String {
        describe_isogeny_frobenius_relation(self)
    }
}

/// Describes the Frobenius data attached to one graph node representative.
pub fn describe_isogeny_graph_node_frobenius_data(node: &IsogenyGraphNodeFrobeniusData) -> String {
    [
        "Isogeny-graph node Frobenius data".to_string(),
        format!("node id: {}", node.node_id().0),
        format!("base field: {}", node.frobenius_trace().base_field()),
        format!("curve order: {}", node.frobenius_trace().curve_order()),
        format!("trace: {}", node.frobenius_trace().trace()),
        format!(
            "same order as reference: {}",
            yes_no(node.same_curve_order_as_reference())
        ),
        format!(
            "same trace as reference: {}",
            yes_no(node.same_trace_as_reference())
        ),
        format!("holds: {}", yes_no(node.holds())),
    ]
    .join("\n")
}

impl Visualizable for IsogenyGraphNodeFrobeniusData {
    fn format_compact(&self) -> String {
        format!("node {}: {}", self.node_id().0, yes_no(self.holds()))
    }

    fn describe(&self) -> String {
        describe_isogeny_graph_node_frobenius_data(self)
    }
}

/// Describes the Frobenius invariance report across an isogeny graph.
pub fn describe_isogeny_graph_frobenius_report(report: &IsogenyGraphFrobeniusReport) -> String {
    let node_lines = report
        .node_reports()
        .iter()
        .map(|node| format!("node {}: {}", node.node_id().0, yes_no(node.holds())))
        .collect::<Vec<_>>()
        .join(", ");

    [
        "Isogeny-graph Frobenius report".to_string(),
        format!("reference node: {}", report.reference_node().0),
        format!("base field: {}", report.reference().base_field()),
        format!("reference order: {}", report.reference().curve_order()),
        format!("reference trace: {}", report.reference().trace()),
        format!("checked nodes: {}", report.checked_nodes()),
        format!("checked edges: {}", report.checked_edges()),
        format!(
            "all same curve order: {}",
            yes_no(report.all_same_curve_order())
        ),
        format!("all same trace: {}", yes_no(report.all_same_trace())),
        format!("holds: {}", yes_no(report.holds())),
        format!("per-node verdicts: {node_lines}"),
    ]
    .join("\n")
}

impl Visualizable for IsogenyGraphFrobeniusReport {
    fn format_compact(&self) -> String {
        format!(
            "graph Frobenius report on {} nodes: {}",
            self.checked_nodes(),
            yes_no(self.holds())
        )
    }

    fn describe(&self) -> String {
        describe_isogeny_graph_frobenius_report(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::frobenius::extension_counts::compare_extension_count_with_enumeration;
    use crate::elliptic_curves::frobenius::group_order::{
        MestreConfig, MestreGroupOrderReport, MestreSide, MestreStepReport,
    };
    use crate::elliptic_curves::frobenius::{
        AbsoluteFrobenius, FrobeniusTrace, RelativeFrobenius,
        characteristic_equation::FrobeniusCharacteristicEquationCurveModel,
        group_order::{GroupOrderReport, SmallFieldGroupOrderStrategy},
        orbit::relative_frobenius_orbit,
    };
    use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
    use crate::elliptic_curves::short_weierstrass::isomorphisms::ShortWeierstrassQuadraticTwist;
    use crate::elliptic_curves::traits::{
        AffineCurveModel, FiniteGroupCurveModel, FrobeniusTraceCurveModel,
    };
    use crate::fields::{Fp, traits::EnumerableFiniteField, traits::Field};
    use crate::isogenies::frobenius_relation::{
        FrobeniusComparableIsogeny, FrobeniusComparableIsogenyGraph,
    };
    use crate::isogenies::graphs::IsogenyGraphBuilder;
    use crate::isogenies::scalar_multiplication::ScalarMultiplicationIsogeny;
    use crate::proptest_support::fields::ProptestF17Sqrt3Field;
    use crate::visualization::elliptic_curves::frobenius::{
        describe_absolute_frobenius, describe_character_sum_point_count,
        describe_frobenius_characteristic_equation_check,
        describe_frobenius_characteristic_equation_exhaustive_report,
        describe_frobenius_characteristic_polynomial, describe_frobenius_extension_count_report,
        describe_frobenius_extension_count_sequence_report,
        describe_frobenius_extension_enumeration_comparison_report,
        describe_frobenius_local_zeta_function, describe_frobenius_on_exact_torsion_report,
        describe_frobenius_orbit, describe_frobenius_trace, describe_group_order_report,
        describe_hasse_bound_report, describe_hasse_interval,
        describe_hasse_multiple_search_report, describe_isogeny_frobenius_relation,
        describe_isogeny_graph_frobenius_report, describe_mestre_group_order_report,
        describe_mestre_step_report, describe_quadratic_twist_frobenius_relation,
        describe_relative_frobenius, format_absolute_frobenius, format_character_sum_point_count,
        format_frobenius_trace, format_group_order_report, format_hasse_interval,
        format_hasse_multiple_search_report, format_mestre_group_order_report,
        format_mestre_step_report, format_relative_frobenius,
    };
    use crate::visualization::traits::Visualizable;
    use num_bigint::BigUint;

    type F17 = Fp<17>;
    type F19 = Fp<19>;
    type F7 = Fp<7>;
    type F41 = Fp<41>;
    type F43 = Fp<43>;
    type F17Squared = ProptestF17Sqrt3Field;

    crate::fields::extension_field::define_fp_quadratic_extension!(
        spec: VisualizationF43Sqrt2Spec,
        field: VisualizationF43Sqrt2,
        base: F43,
        non_residue: 2,
        name: "visualization F43(sqrt(2))",
    );

    fn f41_curve() -> ShortWeierstrassCurve<F41> {
        ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3))
            .expect("valid F41 curve")
    }

    fn lift_f17_curve_to_f17_squared(
        curve: &ShortWeierstrassCurve<F17>,
    ) -> ShortWeierstrassCurve<F17Squared> {
        ShortWeierstrassCurve::<F17Squared>::new(
            F17Squared::from_base(*curve.a()),
            F17Squared::from_base(*curve.b()),
        )
        .expect("lifting a smooth F17 curve should preserve smoothness")
    }

    fn first_nonsquare<F>() -> F::Elem
    where
        F: EnumerableFiniteField + crate::fields::traits::SqrtField,
    {
        F::elements()
            .into_iter()
            .find(|value| !F::is_zero(value) && !F::has_square_root(value))
            .expect("small odd prime fields should contain non-squares")
    }

    fn first_non_fixed_point(
        curve: &ShortWeierstrassCurve<VisualizationF43Sqrt2>,
    ) -> crate::elliptic_curves::AffinePoint<VisualizationF43Sqrt2> {
        for x in VisualizationF43Sqrt2::elements() {
            for y in VisualizationF43Sqrt2::elements() {
                if let Ok(point) = curve.point(x.clone(), y) {
                    let image = curve
                        .absolute_frobenius_power_point(&point, 1)
                        .expect("absolute Frobenius should evaluate");
                    if image != point {
                        return point;
                    }
                }
            }
        }

        panic!("expected a non-fixed point over F43^2")
    }

    #[test]
    fn frobenius_metadata_visualizations_keep_absolute_and_relative_distinct() {
        let absolute = AbsoluteFrobenius::for_field::<F43>(3);
        let relative = RelativeFrobenius::for_field::<F17Squared>(2);

        assert_eq!(format_absolute_frobenius(&absolute), "π_43^3");
        assert_eq!(format_relative_frobenius(&relative), "π_(17^2)^2");

        let absolute_description = describe_absolute_frobenius(&absolute);
        let relative_description = describe_relative_frobenius(&relative);

        assert!(absolute_description.contains("Absolute Frobenius"));
        assert!(absolute_description.contains("characteristic p: 43"));
        assert!(relative_description.contains("Relative Frobenius"));
        assert!(relative_description.contains("base field: F_(17^2)"));
        assert!(relative_description.contains("field order q: 289"));
    }

    #[test]
    fn trace_polynomial_and_zeta_visualizations_share_the_same_frobenius_story() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let trace = curve.frobenius_trace().expect("trace should compute");
        let polynomial = trace.characteristic_polynomial();
        let zeta = trace.local_zeta_function();

        assert!(format_frobenius_trace(&trace).contains("t ="));

        let trace_description = describe_frobenius_trace(&trace);
        let polynomial_description = describe_frobenius_characteristic_polynomial(&polynomial);
        let zeta_description = describe_frobenius_local_zeta_function(&zeta);

        assert!(trace_description.contains("curve order #E(F_q)"));
        assert!(polynomial_description.contains("χ_π(T)"));
        assert!(polynomial_description.contains("discriminant t^2 - 4q"));
        assert!(zeta_description.contains("source characteristic polynomial"));
        assert!(zeta_description.contains("numerator:"));
        assert!(zeta_description.contains("denominator:"));
        assert_eq!(zeta.format_compact(), zeta.pretty());
    }

    #[test]
    fn hasse_and_curve_type_visualizations_explain_their_exact_criteria() {
        let ordinary_curve =
            ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");
        let _ordinary_trace = ordinary_curve
            .frobenius_trace()
            .expect("trace should compute");
        let hasse_report = ordinary_curve
            .verify_hasse_bound()
            .expect("Hasse report should compute");

        let hasse_description = describe_hasse_bound_report(&hasse_report);

        assert!(hasse_description.contains("trace square t^2"));
        assert!(hasse_description.contains("bound square 4q"));
        assert!(hasse_description.contains("slack 4q - t^2"));
    }

    #[test]
    fn hasse_interval_visualization_reports_discrete_search_data() {
        let interval = crate::elliptic_curves::frobenius::HasseInterval::for_q(43)
            .expect("q = 43 should define a Hasse interval");

        assert_eq!(format_hasse_interval(&interval), "H(43) = [31 , 57]");

        let description = describe_hasse_interval(&interval);
        assert!(description.contains("field order q: 43"));
        assert!(description.contains("interval H(q): [31 , 57]"));
        assert!(description.contains("integer candidate count: 27"));
        assert!(description.contains("floor(sqrt(4q)): 13"));
    }

    #[test]
    fn naive_hasse_multiple_search_visualization_reports_the_first_hit_and_steps() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let point = curve
            .generator()
            .expect("the sample curve should be cyclic");
        let report = curve
            .find_annihilating_multiple_in_hasse_interval_naive(&point)
            .expect("naive Hasse search should succeed");

        assert_eq!(
            format_hasse_multiple_search_report(&report),
            "first H(q)-multiple annihilating P: 6"
        );

        let description = describe_hasse_multiple_search_report(&report);
        assert!(description.contains("searched interval: H(7) = [3 , 13]"));
        assert!(description.contains("tested candidates: 4"));
        assert!(description.contains("first annihilating multiple: 6"));
        assert!(description.contains("M = 6 gives [M]P = O"));
    }

    #[test]
    fn character_sum_visualization_reports_the_counting_formula() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let report = curve
            .group_order_by_quadratic_character()
            .expect("character-sum count should succeed");

        assert_eq!(
            format_character_sum_point_count(&report),
            "#E(F_43) via χ-sum = 34"
        );

        let description = describe_character_sum_point_count(&report);
        assert!(description.contains("Quadratic-character point count"));
        assert!(description.contains("character sum Σ χ(f(x))"));
        assert!(description.contains("counting formula: #E(F_q) = q + 1 + Σ χ(f(x))"));
    }

    #[test]
    fn unified_group_order_visualization_mentions_the_strategy() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let report = curve
            .group_order_by_small_field(SmallFieldGroupOrderStrategy::Auto)
            .expect("automatic group order should succeed");

        assert_eq!(
            format_group_order_report(&report),
            "#E(F_43) via χ-sum = 34"
        );

        let description = describe_group_order_report(&report);
        assert!(description.contains("Group order"));
        assert!(description.contains("strategy: quadratic character"));
        assert!(description.contains("curve order #E(F_q): 34"));
    }

    #[test]
    fn mestre_visualizations_show_side_history_and_lower_bounds() {
        let base_field = crate::fields::finite_field_descriptor::FiniteFieldDescriptor::new(
            43,
            core::num::NonZeroU32::new(1).expect("1 is non-zero"),
        )
        .expect("prime field descriptor should build");
        let original = FrobeniusTrace::from_order(base_field.clone(), 52)
            .expect("original Frobenius package should build");
        let twist = FrobeniusTrace::from_order(base_field, 36).expect("twist package should build");
        let point_order_report = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid sample curve")
            .point_order_from_multiple(
                &ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
                    .expect("valid sample curve")
                    .point(F7::from_i64(2), F7::from_i64(1))
                    .expect("sample point should lie on the curve"),
                BigUint::from(6u8),
                &[(BigUint::from(2u8), 1), (BigUint::from(3u8), 1)],
            )
            .expect("known-multiple route should recover a sample order");
        let step = MestreStepReport::new(
            MestreSide::QuadraticTwist,
            45,
            point_order_report,
            BigUint::from(9u8),
        );
        let mestre = MestreGroupOrderReport::new(
            MestreConfig::with_iteration_cap(8),
            MestreSide::QuadraticTwist,
            original,
            twist,
            vec![step.clone()],
        );

        assert_eq!(
            format_mestre_step_report(&step),
            "quadratic twist: M = 45, ord(P) = 6, running λ lower bound = 9"
        );
        assert_eq!(
            format_mestre_group_order_report(&mestre),
            "#E(F_43) via Mestre = 52"
        );

        let step_description = describe_mestre_step_report(&step);
        let mestre_description = describe_mestre_group_order_report(&mestre);
        let unified_description =
            describe_group_order_report(&GroupOrderReport::MestreFp(Box::new(mestre.clone())));

        assert!(step_description.contains("side: quadratic twist"));
        assert!(step_description.contains("annihilating multiple in H(p): 45"));
        assert!(mestre_description.contains("resolved side: quadratic twist"));
        assert!(mestre_description.contains("shared Hasse interval: H(43)"));
        assert!(mestre_description.contains("step 1: quadratic twist: M = 45"));
        assert!(mestre_description.contains("iteration cap: 8"));
        assert!(unified_description.contains("strategy: Mestre"));
        assert!(unified_description.contains("note: the returned group order is always #E(F_p)"));
    }

    #[test]
    fn extension_count_visualizations_show_the_derived_and_exhaustive_routes() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
            .expect("valid F17 curve");
        let lifted_curve = lift_f17_curve_to_f17_squared(&curve);
        let trace = curve.frobenius_trace().expect("trace should compute");
        let report =
            trace.curve_order_over_extension(core::num::NonZeroU32::new(2).expect("2 is positive"));
        let sequence = trace.curve_orders_over_extensions_through(
            core::num::NonZeroU32::new(3).expect("3 is positive"),
        );
        let comparison = compare_extension_count_with_enumeration(&lifted_curve, &trace)
            .expect("comparison should compute");

        let report_description = describe_frobenius_extension_count_report(&report);
        let sequence_description = describe_frobenius_extension_count_sequence_report(&sequence);
        let comparison_description =
            describe_frobenius_extension_enumeration_comparison_report(&comparison);

        assert!(report_description.contains("extension field: F_(17^2)"));
        assert!(report_description.contains("power sum s_n = α^n + β^n"));
        assert!(sequence_description.contains("degree 3"));
        assert!(comparison_description.contains("Frobenius-derived count"));
        assert!(comparison_description.contains("exhaustive count"));
        assert!(comparison_description.contains("agreement: yes"));
    }

    #[test]
    fn characteristic_equation_visualizations_show_all_pointwise_terms() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("sample point should lie on the curve");
        let polynomial = curve
            .frobenius_trace()
            .expect("trace should compute")
            .characteristic_polynomial();
        let check = curve
            .verify_frobenius_characteristic_equation_at_point(&point, &polynomial)
            .expect("pointwise check should compute");
        let exhaustive = curve
            .verify_frobenius_characteristic_equation_exhaustive()
            .expect("exhaustive report should compute");

        let check_description = describe_frobenius_characteristic_equation_check(&check);
        let exhaustive_description =
            describe_frobenius_characteristic_equation_exhaustive_report(&exhaustive);

        assert!(check_description.contains("π_q(P)"));
        assert!(check_description.contains("π_q^2(P)"));
        assert!(check_description.contains("[t]π_q(P)"));
        assert!(check_description.contains("[q]P"));
        assert!(exhaustive_description.contains("checked points"));
        assert!(exhaustive_description.contains("failed checks: 0"));
        assert!(exhaustive_description.contains("failed points: none"));
    }

    #[test]
    fn orbit_and_torsion_visualizations_report_motion_and_periods() {
        let curve = ShortWeierstrassCurve::<VisualizationF43Sqrt2>::new(
            VisualizationF43Sqrt2::from_base(F43::zero()),
            VisualizationF43Sqrt2::from_base(F43::one()),
        )
        .expect("valid extension curve");
        let point = first_non_fixed_point(&curve);
        let orbit = curve
            .absolute_frobenius_orbit(&point, 1)
            .expect("orbit should compute");
        let torsion_report = curve
            .absolute_frobenius_on_exact_torsion(4, 1)
            .expect("torsion report should compute");

        let orbit_description = describe_frobenius_orbit(&orbit);
        let torsion_description = describe_frobenius_on_exact_torsion_report(&torsion_report);

        assert!(orbit_description.contains("period: 2"));
        assert!(orbit_description.contains("points: ["));
        assert!(torsion_description.contains("exact order n: 4"));
        assert!(torsion_description.contains("fixed count: 0"));
        assert!(torsion_description.contains("moved count: 12"));
        assert!(torsion_description.contains("orbit periods: [2, 2, 2, 2, 2, 2]"));
        assert!(
            torsion_description.contains("minimal absolute-Frobenius fixing powers: d = 2: 12")
        );
    }

    #[test]
    fn quadratic_twist_and_isogeny_visualizations_report_their_invariants() {
        let twist_curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))
            .expect("valid curve");
        let twist_package =
            ShortWeierstrassQuadraticTwist::new(twist_curve, first_nonsquare::<F19>())
                .expect("quadratic twist should build");
        let twist_relation = twist_package
            .frobenius_relation()
            .expect("twist relation should compute");

        let isogeny =
            ScalarMultiplicationIsogeny::new(f41_curve(), 2).expect("scalar isogeny should build");
        let isogeny_relation = isogeny
            .frobenius_relation_report()
            .expect("relation should compute");

        let twist_description = describe_quadratic_twist_frobenius_relation(&twist_relation);
        let isogeny_description = describe_isogeny_frobenius_relation(&isogeny_relation);

        assert!(twist_description.contains("expected sum 2q + 2"));
        assert!(twist_description.contains("trace negation t' = -t holds: yes"));
        assert!(isogeny_description.contains("isogeny degree: 4"));
        assert!(isogeny_description.contains("same curve order: yes"));
        assert!(isogeny_description.contains("same trace: yes"));
    }

    #[test]
    fn graph_visualization_reports_reference_and_per_node_verdicts() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("small graph should build");
        let report = graph
            .frobenius_relation_report()
            .expect("graph report should compute");

        let description = describe_isogeny_graph_frobenius_report(&report);

        assert!(description.contains("reference node: 0"));
        assert!(description.contains("checked nodes:"));
        assert!(description.contains("checked edges:"));
        assert!(description.contains("per-node verdicts:"));
        assert!(description.contains("node 0: yes"));
    }

    #[test]
    fn visualizable_trait_is_hooked_up_for_frobenius_objects() {
        let relative = RelativeFrobenius::for_field::<F17Squared>(1);
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("sample point should lie on the curve");
        let orbit =
            relative_frobenius_orbit(&curve, &point).expect("relative orbit should compute");

        assert!(relative.format_compact().contains("π_(17^2)"));
        assert!(relative.describe().contains("field order q: 289"));
        assert!(orbit.format_compact().contains("period 1 orbit"));
        assert!(orbit.describe().contains("Frobenius orbit"));
    }
}
