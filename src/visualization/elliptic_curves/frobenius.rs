use std::collections::BTreeMap;

use num_bigint::BigUint;

use crate::elliptic_curves::frobenius::{
    AbsoluteFrobenius, FrobeniusCharacteristicPolynomial, FrobeniusCurveType,
    FrobeniusLocalZetaFunction, FrobeniusTrace, HasseInterval, RelativeFrobenius,
    character_sum::CharacterSumPointCount,
    characteristic_equation::{
        FrobeniusCharacteristicEquationCheck, FrobeniusCharacteristicEquationExhaustiveReport,
    },
    extension_counts::{
        FrobeniusExtensionCountReport, FrobeniusExtensionCountSequenceReport,
        FrobeniusExtensionEnumerationComparisonReport,
    },
    group_order::{
        GroupOrderReport, GroupOrderRoute, MestreGroupOrderReport, MestreSide, MestreStepReport,
    },
    hasse::{HasseBoundReport, HasseMultipleSearchReport, HasseMultipleSearchStep},
    orbit::FrobeniusOrbit,
    quadratic_twist::QuadraticTwistFrobeniusRelation,
    schoof::{
        SchoofGroupOrderOutcome, SchoofGroupOrderReport, SchoofTraceCrtOutcome,
        SchoofTraceCrtReport, SchoofTraceModOddPrimeOutcome,
    },
    torsion::{FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport},
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::frobenius_relation::{
    IsogenyFrobeniusRelation, IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData,
};
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::short_weierstrass::{
        describe_point_order_from_multiple_report, format_point_order_from_multiple_report,
    },
    polynomials::dense::format_dense_polynomial,
    shared::{comma_list, compact_visualizable_list, yes_no},
};

fn field_symbol(characteristic: &BigUint, extension_degree: u32) -> String {
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
        format!("[{}]", compact_visualizable_list(points.iter()))
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

    Some(comma_list(
        histogram
            .into_iter()
            .map(|(power, count)| format!("d = {power}: {count}")),
    ))
}

/// Formats the absolute Frobenius metadata compactly.
fn format_absolute_frobenius(frobenius: &AbsoluteFrobenius) -> String {
    iterated_symbol(
        &format!("π_{}", frobenius.characteristic()),
        frobenius.power(),
    )
}

/// Describes the absolute Frobenius metadata.
fn describe_absolute_frobenius(frobenius: &AbsoluteFrobenius) -> String {
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
fn format_relative_frobenius(frobenius: &RelativeFrobenius) -> String {
    let symbol = field_symbol(
        &frobenius.base_field().characteristic,
        frobenius.base_field().extension_degree.get(),
    );
    iterated_symbol(&format!("π_{symbol}"), frobenius.power())
}

/// Describes the relative Frobenius metadata.
fn describe_relative_frobenius(frobenius: &RelativeFrobenius) -> String {
    [
        "Relative Frobenius".to_string(),
        format!("symbol: {}", format_relative_frobenius(frobenius)),
        format!("base field: {}", frobenius.base_field()),
        format!(
            "field order q: {}",
            frobenius.base_field().cardinality()
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
fn format_frobenius_trace(trace: &FrobeniusTrace) -> String {
    format!("t = {} over {}", trace.trace(), trace.base_field())
}

/// Describes a Frobenius trace package recovered from `#E(F_q)`.
fn describe_frobenius_trace(trace: &FrobeniusTrace) -> String {
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
pub(crate) fn format_hasse_interval(interval: &HasseInterval) -> String {
    format!(
        "H({}) = [{} , {}]",
        interval.q(),
        interval.lower(),
        interval.upper()
    )
}

/// Describes the discrete Hasse interval of possible values of `#E(F_q)`.
fn describe_hasse_interval(interval: &HasseInterval) -> String {
    let doubled_sqrt_floor = interval.upper() - (interval.q() + BigUint::from(1u8));
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
fn format_hasse_multiple_search_step<P: Visualizable>(step: &HasseMultipleSearchStep<P>) -> String {
    format!(
        "M = {} gives [M]P = {}",
        step.candidate_multiple(),
        step.image().format_compact()
    )
}

/// Describes one tested candidate in a Hasse-interval multiple search.
fn describe_hasse_multiple_search_step<P: Visualizable>(
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
pub(crate) fn format_hasse_multiple_search_report<P: Visualizable>(
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
pub(crate) fn describe_hasse_multiple_search_report<P: Visualizable>(
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

mod group_order;

use group_order::{
    describe_character_sum_point_count, format_character_sum_point_count, format_group_order_report,
};
pub(crate) use group_order::{
    describe_group_order_report, describe_mestre_step_report, format_mestre_step_report,
};

#[cfg(test)]
use group_order::{describe_mestre_group_order_report, format_mestre_group_order_report};

/// Describes the Frobenius characteristic polynomial.
fn describe_frobenius_characteristic_polynomial(
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
fn describe_frobenius_local_zeta_function(zeta: &FrobeniusLocalZetaFunction) -> String {
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
fn describe_hasse_bound_report(report: &HasseBoundReport) -> String {
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
fn format_frobenius_orbit<P: Visualizable>(orbit: &FrobeniusOrbit<P>) -> String {
    format!(
        "period {} orbit {}",
        orbit.period(),
        orbit_points_text(orbit.points())
    )
}

/// Describes one Frobenius orbit.
fn describe_frobenius_orbit<P: Visualizable>(orbit: &FrobeniusOrbit<P>) -> String {
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
        compact_visualizable_list(report.points().iter())
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
fn describe_frobenius_extension_count_report(report: &FrobeniusExtensionCountReport) -> String {
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
fn describe_frobenius_extension_count_sequence_report(
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
fn describe_frobenius_extension_enumeration_comparison_report(
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
fn describe_frobenius_characteristic_equation_check<P: Visualizable>(
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
fn describe_frobenius_characteristic_equation_exhaustive_report<P: Visualizable>(
    report: &FrobeniusCharacteristicEquationExhaustiveReport<P>,
) -> String {
    let failed_points = if report.failed_checks().is_empty() {
        "none".to_string()
    } else {
        comma_list(
            report
                .failed_checks()
                .iter()
                .map(|check| check.point().format_compact()),
        )
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
fn describe_quadratic_twist_frobenius_relation(
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
fn describe_isogeny_frobenius_relation(relation: &IsogenyFrobeniusRelation) -> String {
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
fn describe_isogeny_graph_node_frobenius_data(node: &IsogenyGraphNodeFrobeniusData) -> String {
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
fn describe_isogeny_graph_frobenius_report(report: &IsogenyGraphFrobeniusReport) -> String {
    let node_lines = comma_list(
        report
            .node_reports()
            .iter()
            .map(|node| format!("node {}: {}", node.node_id().0, yes_no(node.holds()))),
    );

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
mod tests;
