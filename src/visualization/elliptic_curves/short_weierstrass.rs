use core::fmt;

use num_bigint::BigUint;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::short_weierstrass::{
    ExhaustivePointOrderReport, ExponentAccumulationReport, ExponentAccumulationStep,
    ExponentLowerBoundGroupOrderVerification, GroupExponentReport, GroupExponentStrategy,
    HasseIntervalPointOrderReport, PointOrderFromMultipleReport, PointOrderReductionStep,
    PointOrderReport, PointOrderStrategyKind, ShortWeierstrassCurve,
};
use crate::elliptic_curves::traits::{
    CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure, FiniteGroupCurveModel,
    GroupCurveModel,
};
use crate::fields::{EnumerableFiniteField, Field, SqrtField};
use crate::visualization::elliptic_curves::frobenius::{
    describe_group_order_report, describe_hasse_multiple_search_report, format_hasse_interval,
    format_hasse_multiple_search_report,
};
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::traits::Visualizable;

fn format_elem<F>(value: &F::Elem) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    value.format_elem()
}

fn parenthesize_if_needed(text: &str) -> String {
    if text.contains(' ') || text.starts_with('-') || text.contains('/') {
        format!("({text})")
    } else {
        text.to_string()
    }
}

fn equation_string<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let mut terms = vec!["x^3".to_string()];

    if !F::is_zero(curve.a()) {
        if F::eq(curve.a(), &F::one()) {
            terms.push("x".to_string());
        } else {
            terms.push(format!(
                "{}x",
                parenthesize_if_needed(&format_elem::<F>(curve.a()))
            ));
        }
    }

    if !F::is_zero(curve.b()) {
        terms.push(parenthesize_if_needed(&format_elem::<F>(curve.b())));
    }

    format!("y^2 = {}", terms.join(" + "))
}

/// Formats a short-Weierstrass curve compactly.
pub fn format_curve<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    equation_string(curve)
}

/// Formats an affine point compactly.
pub fn format_point<F>(point: &AffinePoint<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    point.to_coordinates_string()
}

/// Formats an affine point using the compact field-element visualization.
pub fn format_point_compact<F>(point: &AffinePoint<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    match point {
        AffinePoint::Infinity => "O".to_string(),
        AffinePoint::Finite { x, y } => {
            format!("({}, {})", format_elem::<F>(x), format_elem::<F>(y))
        }
    }
}

impl<F> Visualizable for AffinePoint<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        match self {
            AffinePoint::Infinity => "O".to_string(),
            AffinePoint::Finite { x, y } => {
                format!("({}, {})", format_elem::<F>(x), format_elem::<F>(y))
            }
        }
    }

    fn describe(&self) -> String {
        match self {
            AffinePoint::Infinity => {
                "Affine point\npoint: O\nrole: distinguished identity point".to_string()
            }
            AffinePoint::Finite { x, y } => format!(
                "Affine point\npoint: {}\nx-coordinate: {}\ny-coordinate: {}",
                self.format_compact(),
                format_elem::<F>(x),
                format_elem::<F>(y)
            ),
        }
    }
}

impl<F> Visualizable for ShortWeierstrassCurve<F>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    fn format_compact(&self) -> String {
        format_curve(self)
    }

    fn describe(&self) -> String {
        describe_curve(self)
    }
}

/// Describes a short-Weierstrass curve with its standard invariants.
pub fn describe_curve<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    [
        "Short-Weierstrass curve".to_string(),
        format!("equation: {}", equation_string(curve)),
        format!("a: {}", format_elem::<F>(curve.a())),
        format!("b: {}", format_elem::<F>(curve.b())),
        format!("discriminant: {}", format_elem::<F>(&curve.discriminant())),
        format!("c4: {}", format_elem::<F>(&curve.c4())),
        format!("c6: {}", format_elem::<F>(&curve.c6())),
        format!("j-invariant: {}", format_elem::<F>(&curve.j_invariant())),
    ]
    .join("\n")
}

/// Describes a point together with its role relative to a chosen curve.
pub fn describe_point<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let mut lines = vec![
        "Curve point".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("point: {}", format_point(point)),
        format!(
            "identity: {}",
            if curve.is_identity(point) {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "on curve: {}",
            if curve.contains(point) { "yes" } else { "no" }
        ),
    ];

    if let AffinePoint::Finite { x, y } = point {
        lines.push(format!("x-coordinate: {}", format_elem::<F>(x)));
        lines.push(format!("y-coordinate: {}", format_elem::<F>(y)));
    }

    lines.join("\n")
}

/// Explains curve membership by comparing `y^2` with `x^3 + ax + b`.
pub fn describe_membership<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    match point {
        AffinePoint::Infinity => [
            "Curve membership".to_string(),
            format!("curve: {}", equation_string(curve)),
            "point: O".to_string(),
            "convention: the point at infinity is part of the curve model".to_string(),
            "result: on curve".to_string(),
        ]
        .join("\n"),
        AffinePoint::Finite { x, y } => {
            let left = F::square(y);
            let x_cubed = F::cube(x);
            let ax = F::mul(curve.a(), x);
            let right = F::add(&F::add(&x_cubed, &ax), curve.b());
            let verdict = if F::eq(&left, &right) {
                "on curve"
            } else {
                "not on curve"
            };

            [
                "Curve membership".to_string(),
                format!("curve: {}", equation_string(curve)),
                format!("point: {}", format_point(point)),
                format!(
                    "left side: y^2 = {}^2 = {}",
                    format_elem::<F>(y),
                    format_elem::<F>(&left)
                ),
                format!(
                    "right side: x^3 + ax + b = {} + {} + {} = {}",
                    format_elem::<F>(&x_cubed),
                    format_elem::<F>(&ax),
                    format_elem::<F>(curve.b()),
                    format_elem::<F>(&right)
                ),
                format!(
                    "comparison: {} = {}",
                    format_elem::<F>(&left),
                    format_elem::<F>(&right)
                ),
                format!("result: {verdict}"),
            ]
            .join("\n")
        }
    }
}

/// Explains affine point addition on a short-Weierstrass curve.
pub fn explain_add<F>(
    curve: &ShortWeierstrassCurve<F>,
    left: &AffinePoint<F>,
    right: &AffinePoint<F>,
) -> Result<String, CurveError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    if !curve.contains(left) || !curve.contains(right) {
        return Err(CurveError::PointNotOnCurve);
    }

    let mut lines = vec![
        "Point addition".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("left: {}", format_point(left)),
        format!("right: {}", format_point(right)),
    ];

    let result = match (left, right) {
        (AffinePoint::Infinity, _) => {
            lines.push("case: O + Q = Q".to_string());
            right.clone()
        }
        (_, AffinePoint::Infinity) => {
            lines.push("case: P + O = P".to_string());
            left.clone()
        }
        (
            AffinePoint::Finite {
                x: x_left,
                y: y_left,
            },
            AffinePoint::Finite {
                x: x_right,
                y: y_right,
            },
        ) => {
            if F::eq(x_left, x_right) && F::is_zero(&F::add(y_left, y_right)) {
                lines.push("case: inverse points with the same x-coordinate".to_string());
                lines.push("result: O".to_string());
                curve.identity()
            } else if F::eq(x_left, x_right) {
                let numerator = F::add(&F::mul(&F::from_i64(3), &F::square(x_left)), curve.a());
                let denominator = F::mul(&F::from_i64(2), y_left);
                let slope = F::div(&numerator, &denominator)
                    .expect("doubling denominator is non-zero in this branch");
                let doubled = curve.double(left)?;

                lines.push("case: tangent formula for doubling".to_string());
                lines.push(format!(
                    "slope: (3x^2 + a) / (2y) = {} / {} = {}",
                    format_elem::<F>(&numerator),
                    format_elem::<F>(&denominator),
                    format_elem::<F>(&slope)
                ));
                doubled
            } else {
                let numerator = F::sub(y_right, y_left);
                let denominator = F::sub(x_right, x_left);
                let slope = F::div(&numerator, &denominator)
                    .expect("addition denominator is non-zero in this branch");
                let sum = curve.add(left, right)?;

                lines.push("case: secant formula for distinct x-coordinates".to_string());
                lines.push(format!(
                    "slope: (y2 - y1) / (x2 - x1) = {} / {} = {}",
                    format_elem::<F>(&numerator),
                    format_elem::<F>(&denominator),
                    format_elem::<F>(&slope)
                ));
                sum
            }
        }
    };

    lines.push(format!("result: {}", format_point(&result)));
    Ok(lines.join("\n"))
}

/// Lists every point of a small finite curve group.
pub fn list_points<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let points = curve.points();
    let mut lines = vec![
        "Curve points".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("group order: {}", points.len()),
    ];

    for (index, point) in points.iter().enumerate() {
        lines.push(format!("{index}: {}", format_point(point)));
    }

    lines.join("\n")
}

/// Describes the order of a point in a small finite curve group.
pub fn describe_point_order<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let mut lines = vec![
        "Point order".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("point: {}", format_point(point)),
        "method: repeated addition up to the full enumerated group order".to_string(),
    ];

    match curve.point_order(point) {
        Some(order) => {
            lines.push(format!("group order: {}", curve.order()));
            lines.push(format!("point order: {order}"));
        }
        None => {
            lines.push("result: point is not on the curve".to_string());
        }
    }

    lines.join("\n")
}

fn format_invariant_factor_surface(structure: FiniteAbelianGroupStructure) -> String {
    if structure.order == 1 {
        return "trivial group".to_string();
    }

    if structure.cyclic {
        return format!("Z/{}Z", structure.order);
    }

    match structure.invariant_factors {
        Some((left, right)) => format!("Z/{left}Z x Z/{right}Z"),
        None => format!("order {}, exponent {}", structure.order, structure.exponent),
    }
}

/// Describes the finite abelian group structure of a small enumerated curve.
pub fn describe_group_structure<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let structure = curve.group_structure();

    [
        "Finite curve group structure".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("group order: {}", structure.order),
        format!("cyclic: {}", if structure.cyclic { "yes" } else { "no" }),
        format!("exponent: {}", structure.exponent),
        format!(
            "invariant factors: {}",
            format_invariant_factor_surface(structure)
        ),
    ]
    .join("\n")
}

/// Returns a compact educational summary of the finite group structure.
pub fn summarize_group_structure<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let structure = curve.group_structure();

    [
        format!("cyclic: {}", if structure.cyclic { "yes" } else { "no" }),
        format!("exponent: {}", structure.exponent),
        format!(
            "invariant factors: {}",
            format_invariant_factor_surface(structure)
        ),
    ]
    .join("\n")
}

/// Describes how many points have each exact order on a small finite curve.
pub fn describe_order_distribution<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let distribution = curve.order_distribution();
    let mut lines = vec![
        "Point-order distribution".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("group order: {}", curve.order()),
    ];

    if distribution.is_empty() {
        lines.push("distribution: no enumerated points".to_string());
        return lines.join("\n");
    }

    for (order, count) in distribution {
        lines.push(format!("order {order}: {count} point(s)"));
    }

    lines.join("\n")
}

/// Returns a compact point-order distribution summary.
pub fn summarize_order_distribution<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    curve
        .order_distribution()
        .into_iter()
        .map(|(order, count)| format!("{order} -> {count}"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Explains scalar multiplication on a curve point.
pub fn describe_scalar_mul<F>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    scalar: i64,
) -> Result<String, CurveError>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let result = curve.mul_scalar_signed(point, scalar)?;
    let magnitude = scalar.unsigned_abs();
    let mut lines = vec![
        "Scalar multiplication".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("point: {}", format_point(point)),
        format!("scalar: {scalar}"),
    ];

    if scalar == 0 {
        lines.push("case: [0]P is the identity by definition".to_string());
    } else if scalar < 0 {
        lines.push(format!(
            "case: [{}]P = [{}](-P) with -P = {}",
            scalar,
            magnitude,
            format_point(&curve.neg(point))
        ));
    } else {
        lines.push(format!(
            "method: double-and-add for the binary expansion of {}",
            scalar
        ));
    }

    lines.push(format!("result: [{}]P = {}", scalar, format_point(&result)));
    Ok(lines.join("\n"))
}

/// Explains why a point has its exact order in a small finite curve group.
pub fn explain_point_order<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let mut lines = vec![
        "Point-order explanation".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("point: {}", format_point(point)),
    ];

    let Some(order) = curve.point_order(point) else {
        lines.push("result: point is not on the curve".to_string());
        return lines.join("\n");
    };

    lines.push(format!("group order: {}", curve.order()));
    lines.push(
        "search: enumerate [n]P for 1 <= n <= #E(F_q) until the identity appears".to_string(),
    );

    let mut multiple = curve.identity();
    for n in 1..=order {
        multiple = curve
            .add(&multiple, point)
            .expect("on-curve point should add successfully during order explanation");
        lines.push(format!("[{n}]P = {}", format_point(&multiple)));
    }

    lines.push(format!("first identity hit: [{}]P = O", order));
    lines.push(format!("point order: {order}"));
    lines.join("\n")
}

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
pub fn format_point_order_from_multiple_report(report: &PointOrderFromMultipleReport) -> String {
    format!(
        "ord(P) from M = {} is {}",
        report.supplied_multiple(),
        report.exact_order()
    )
}

/// Describes the prime-peeling order recovery from one known multiple.
pub fn describe_point_order_from_multiple_report(report: &PointOrderFromMultipleReport) -> String {
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

fn point_order_strategy_kind_label(strategy: PointOrderStrategyKind) -> &'static str {
    match strategy {
        PointOrderStrategyKind::Exhaustive => "exhaustive",
        PointOrderStrategyKind::FromKnownMultiple => "from known multiple",
        PointOrderStrategyKind::HasseIntervalNaive => "naive Hasse interval",
    }
}

fn group_exponent_strategy_label(strategy: &GroupExponentStrategy) -> &'static str {
    match strategy {
        GroupExponentStrategy::Exhaustive => "exhaustive",
        GroupExponentStrategy::RandomPoints { .. } => "random points",
    }
}

fn group_order_strategy_label_for_order_route(
    strategy: crate::elliptic_curves::GroupOrderStrategy,
) -> &'static str {
    match strategy {
        crate::elliptic_curves::GroupOrderStrategy::Auto => "auto",
        crate::elliptic_curves::GroupOrderStrategy::Exhaustive => "exhaustive",
        crate::elliptic_curves::GroupOrderStrategy::QuadraticCharacter => "quadratic character",
        crate::elliptic_curves::GroupOrderStrategy::MestreFp(_) => "Mestre",
        crate::elliptic_curves::GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
            ..
        } => "from exponent lower bound",
    }
}

/// Formats an exhaustive point-order report compactly.
pub fn format_exhaustive_point_order_report(report: &ExhaustivePointOrderReport) -> String {
    format!("ord(P) via exhaustive search = {}", report.exact_order())
}

/// Describes an exhaustive point-order report.
pub fn describe_exhaustive_point_order_report(report: &ExhaustivePointOrderReport) -> String {
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
pub fn format_hasse_interval_point_order_report<P: Visualizable>(
    report: &HasseIntervalPointOrderReport<P>,
) -> String {
    format!(
        "ord(P) via H(q) search = {}",
        report.order_from_multiple().exact_order()
    )
}

/// Describes a Hasse-interval point-order report.
pub fn describe_hasse_interval_point_order_report<P: Visualizable>(
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
            group_order_strategy_label_for_order_route(report.group_order_report().strategy())
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
pub fn format_point_order_report<P: Visualizable>(report: &PointOrderReport<P>) -> String {
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
pub fn describe_point_order_report<P: Visualizable>(report: &PointOrderReport<P>) -> String {
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

/// Formats one random-point exponent-accumulation step compactly.
pub fn format_exponent_accumulation_step<P: Visualizable>(
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
pub fn describe_exponent_accumulation_step<P: Visualizable>(
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
pub fn format_exhaustive_group_exponent_report(exact_exponent: &BigUint) -> String {
    format!("group exponent = {exact_exponent}")
}

/// Describes the exact exhaustive group-exponent report.
pub fn describe_exhaustive_group_exponent_report(exact_exponent: &BigUint) -> String {
    [
        "Exhaustive group exponent".to_string(),
        format!("exact exponent: {exact_exponent}"),
        "strategy: compute every point order in the tiny ambient group and take their lcm"
            .to_string(),
    ]
    .join("\n")
}

/// Formats the random-point exponent accumulation report compactly.
pub fn format_exponent_accumulation_report<P: Visualizable>(
    report: &ExponentAccumulationReport<P>,
) -> String {
    format!(
        "group exponent lower bound after {} sample(s) = {}",
        report.samples_taken(),
        report.exponent_lower_bound()
    )
}

/// Describes the random-point exponent accumulation report.
pub fn describe_exponent_accumulation_report<P: Visualizable>(
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
            if report.completed_requested_samples() {
                "yes"
            } else {
                "no"
            }
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
pub fn format_exponent_lower_bound_group_order_verification(
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
pub fn describe_exponent_lower_bound_group_order_verification(
    report: &ExponentLowerBoundGroupOrderVerification,
) -> String {
    let mut lines = vec![
        "Exponent lower-bound verification by group order".to_string(),
        format!("exponent lower bound: {}", report.exponent_lower_bound()),
        format!(
            "group-order route: {}",
            group_order_strategy_label_for_order_route(report.group_order_report().strategy())
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
pub fn format_group_exponent_report<P: Visualizable>(report: &GroupExponentReport<P>) -> String {
    match report {
        GroupExponentReport::Exhaustive(exact_exponent) => {
            format_exhaustive_group_exponent_report(exact_exponent)
        }
        GroupExponentReport::RandomPoints(report) => format_exponent_accumulation_report(report),
    }
}

/// Describes a unified group-exponent report.
pub fn describe_group_exponent_report<P: Visualizable>(report: &GroupExponentReport<P>) -> String {
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

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_bigint::BigUint;
    use num_rational::BigRational;

    use crate::elliptic_curves::{AffineCurveModel, AffinePoint, EnumerableCurveModel};
    use crate::fields::{Field, Fp, Q};
    use crate::visualization::Visualizable;

    use crate::visualization::elliptic_curves::{
        describe_curve, describe_exhaustive_group_exponent_report,
        describe_exhaustive_point_order_report,
        describe_exponent_lower_bound_group_order_verification, describe_group_exponent_report,
        describe_group_structure, describe_membership, describe_order_distribution, describe_point,
        describe_point_order, describe_point_order_from_multiple_report,
        describe_point_order_report, describe_scalar_mul, explain_add, explain_point_order,
        format_curve, format_exponent_lower_bound_group_order_verification,
        format_group_exponent_report, format_point, format_point_compact,
        format_point_order_from_multiple_report, format_point_order_report, list_points,
        summarize_group_structure, summarize_order_distribution,
    };

    type F7 = Fp<7>;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    fn f7_curve() -> crate::elliptic_curves::ShortWeierstrassCurve<F7> {
        crate::elliptic_curves::ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve")
    }

    fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
        f7_curve()
            .point(F7::from_i64(x), F7::from_i64(y))
            .expect("point should lie on the curve")
    }

    #[test]
    fn curve_display_and_equation_string_share_one_equation_surface() {
        let curve = f7_curve();

        assert_eq!(
            curve.to_equation_string(),
            "y^2 = x^3 + (2 (mod 7))x + (3 (mod 7))"
        );
        assert_eq!(format!("{curve}"), curve.to_equation_string());
        assert_eq!(format_curve(&curve), "y^2 = x^3 + 2x + 3");
    }

    #[test]
    fn point_display_uses_affine_coordinates_or_identity_symbol() {
        let point = f7_point(2, 1);
        let infinity = AffinePoint::<F7>::infinity();

        assert_eq!(point.to_coordinates_string(), "(2 (mod 7), 1 (mod 7))");
        assert_eq!(format!("{point}"), point.to_coordinates_string());
        assert_eq!(format_point(&point), point.to_coordinates_string());
        assert_eq!(format_point_compact(&point), "(2, 1)");
        assert_eq!(format_point(&infinity), "O");
        assert_eq!(format_point_compact(&infinity), "O");
    }

    #[test]
    fn debug_output_is_more_informative_than_the_default_derives() {
        let curve = f7_curve();
        let point = f7_point(2, 1);

        assert!(format!("{curve:?}").contains("ShortWeierstrassCurve"));
        assert!(format!("{curve:?}").contains("equation"));
        assert!(format!("{point:?}").contains("AffinePoint"));
        assert!(format!("{point:?}").contains("x"));
    }

    #[test]
    fn curve_description_mentions_invariants() {
        let description = describe_curve(&f7_curve());

        assert!(description.contains("Short-Weierstrass curve"));
        assert!(description.contains("discriminant"));
        assert!(description.contains("j-invariant"));
    }

    #[test]
    fn point_description_mentions_identity_and_membership_status() {
        let description = describe_point(&f7_curve(), &f7_point(2, 1));

        assert!(description.contains("Curve point"));
        assert!(description.contains("identity: no"));
        assert!(description.contains("on curve: yes"));
    }

    #[test]
    fn membership_description_shows_both_sides_of_the_equation() {
        let description = describe_membership(&f7_curve(), &f7_point(2, 1));

        assert!(description.contains("left side: y^2"));
        assert!(description.contains("right side: x^3 + ax + b"));
        assert!(description.contains("result: on curve"));
    }

    #[test]
    fn membership_description_is_honest_about_the_point_at_infinity() {
        let description = describe_membership(&f7_curve(), &AffinePoint::<F7>::infinity());

        assert!(description.contains("point: O"));
        assert!(description.contains("convention"));
    }

    #[test]
    fn addition_explanation_mentions_the_geometric_case_and_result() {
        let explanation =
            explain_add(&f7_curve(), &f7_point(2, 1), &f7_point(3, 1)).expect("valid addition");

        assert!(explanation.contains("Point addition"));
        assert!(explanation.contains("case: secant formula"));
        assert!(explanation.contains("result: (2 (mod 7), 6 (mod 7))"));
    }

    #[test]
    fn point_listing_shows_group_order_and_identity() {
        let listing = list_points(&f7_curve());

        assert!(listing.contains("Curve points"));
        assert!(listing.contains("group order: 6"));
        assert!(listing.contains("0: O"));
    }

    #[test]
    fn point_order_description_mentions_repeated_addition_method() {
        let description = describe_point_order(&f7_curve(), &f7_point(2, 1));

        assert!(description.contains("Point order"));
        assert!(description.contains("repeated addition"));
        assert!(description.contains("point order: 6"));
    }

    #[test]
    fn point_order_description_is_honest_about_invalid_points() {
        let description = describe_point_order(
            &f7_curve(),
            &AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2)),
        );

        assert!(description.contains("result: point is not on the curve"));
    }

    #[test]
    fn group_structure_description_reports_small_cyclic_example() {
        let description = describe_group_structure(&f7_curve());

        assert!(description.contains("Finite curve group structure"));
        assert!(description.contains("group order: 6"));
        assert!(description.contains("cyclic: yes"));
        assert!(description.contains("exponent: 6"));
        assert!(description.contains("invariant factors: Z/6Z"));
    }

    #[test]
    fn compact_group_structure_summary_reports_core_invariants() {
        let summary = summarize_group_structure(&f7_curve());

        assert!(summary.contains("cyclic: yes"));
        assert!(summary.contains("exponent: 6"));
        assert!(summary.contains("invariant factors: Z/6Z"));
    }

    #[test]
    fn order_distribution_description_lists_exact_point_orders() {
        let description = describe_order_distribution(&f7_curve());

        assert!(description.contains("Point-order distribution"));
        assert!(description.contains("order 1: 1 point(s)"));
        assert!(description.contains("order 2: 1 point(s)"));
        assert!(description.contains("order 3: 2 point(s)"));
        assert!(description.contains("order 6: 2 point(s)"));
    }

    #[test]
    fn compact_order_distribution_summary_uses_arrow_surface() {
        let summary = summarize_order_distribution(&f7_curve());

        assert!(summary.contains("1 -> 1"));
        assert!(summary.contains("2 -> 1"));
        assert!(summary.contains("3 -> 2"));
        assert!(summary.contains("6 -> 2"));
    }

    #[test]
    fn scalar_multiplication_description_reports_method_and_result() {
        let description =
            describe_scalar_mul(&f7_curve(), &f7_point(2, 1), 3).expect("valid scalar multiply");

        assert!(description.contains("Scalar multiplication"));
        assert!(description.contains("scalar: 3"));
        assert!(description.contains("double-and-add"));
        assert!(description.contains("result: [3]P = (6 (mod 7), 0 (mod 7))"));
    }

    #[test]
    fn point_order_explanation_lists_successive_multiples_until_identity() {
        let description = explain_point_order(&f7_curve(), &f7_point(2, 1));

        assert!(description.contains("Point-order explanation"));
        assert!(description.contains("[1]P = (2 (mod 7), 1 (mod 7))"));
        assert!(description.contains("[6]P = O"));
        assert!(description.contains("first identity hit: [6]P = O"));
        assert!(description.contains("point order: 6"));
    }

    #[test]
    fn point_order_from_multiple_visualization_reports_the_prime_peeling_steps() {
        let report = f7_curve()
            .point_order_from_multiple(&f7_point(6, 0), bu(6), &[(bu(2), 1), (bu(3), 1)])
            .expect("valid reduction report should build");

        assert_eq!(
            format_point_order_from_multiple_report(&report),
            "ord(P) from M = 6 is 2"
        );

        let description = describe_point_order_from_multiple_report(&report);
        assert!(description.contains("Point order from multiple"));
        assert!(description.contains("supplied multiple M: 6"));
        assert!(description.contains("exact order recovered: 2"));
        assert!(
            description.contains(
                "prime 3: exponent in M = 1, removed exponent = 1, remaining multiple = 2"
            )
        );
    }

    #[test]
    fn unified_point_order_visualization_mentions_the_selected_strategy() {
        let report = f7_curve()
            .point_order_by(
                &f7_point(2, 1),
                crate::elliptic_curves::PointOrderStrategy::HasseIntervalNaive {
                    group_order_strategy: crate::elliptic_curves::GroupOrderStrategy::Auto,
                },
            )
            .expect("Hasse-interval order recovery should succeed");

        assert_eq!(
            format_point_order_report(&report),
            "ord(P) via H(q) search = 6"
        );

        let description = describe_point_order_report(&report);
        assert!(description.contains("Point order report"));
        assert!(description.contains("strategy: naive Hasse interval"));
        assert!(description.contains("exact order: 6"));
        assert!(description.contains("group-order route: quadratic character"));
        assert!(description.contains("first H(q)-multiple annihilating P: 6"));
    }

    #[test]
    fn exhaustive_point_order_visualization_stays_honest_about_the_route() {
        let report = f7_curve()
            .point_order_by(
                &f7_point(2, 1),
                crate::elliptic_curves::PointOrderStrategy::Exhaustive,
            )
            .expect("exhaustive order recovery should succeed");

        let crate::elliptic_curves::PointOrderReport::Exhaustive(exhaustive) = report else {
            panic!("expected the exhaustive route to preserve its variant");
        };

        assert_eq!(
            describe_exhaustive_point_order_report(&exhaustive),
            exhaustive.describe()
        );
        assert!(exhaustive.describe().contains("Exhaustive point order"));
        assert!(exhaustive.describe().contains("exact order: 6"));
    }

    #[test]
    fn group_exponent_visualization_mentions_the_selected_strategy() {
        let curve = f7_curve();
        let sampled_point = f7_point(2, 1);
        let point_index = curve
            .points()
            .iter()
            .position(|candidate| candidate == &sampled_point)
            .expect("sample point should appear in the enumerated group");
        let mut sampler =
            move |upper_bound: usize| (point_index < upper_bound).then_some(point_index);

        let report = curve
            .group_exponent_by(
                crate::elliptic_curves::GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy:
                        crate::elliptic_curves::PointOrderStrategy::HasseIntervalNaive {
                            group_order_strategy: crate::elliptic_curves::GroupOrderStrategy::Auto,
                        },
                },
                &mut sampler,
            )
            .expect("random-point exponent accumulation should succeed");

        assert_eq!(
            format_group_exponent_report(&report),
            "group exponent lower bound after 1 sample(s) = 6"
        );

        let description = describe_group_exponent_report(&report);
        assert!(description.contains("Group exponent report"));
        assert!(description.contains("strategy: random points"));
        assert!(description.contains("exponent lower bound: 6"));
        assert!(description.contains("exact exponent: not certified by this route"));
        assert!(description.contains("point-order route: naive Hasse interval"));
    }

    #[test]
    fn exhaustive_group_exponent_visualization_stays_honest_about_exactness() {
        let curve = f7_curve();
        let mut sampler = |_| Some(0usize);
        let report = curve
            .group_exponent_by(
                crate::elliptic_curves::GroupExponentStrategy::Exhaustive,
                &mut sampler,
            )
            .expect("exhaustive exponent route should succeed");

        let crate::elliptic_curves::GroupExponentReport::Exhaustive(exact_exponent) = report else {
            panic!("expected the exhaustive group-exponent route to preserve its variant");
        };

        assert_eq!(
            describe_exhaustive_group_exponent_report(&exact_exponent),
            "Exhaustive group exponent\nexact exponent: 6\nstrategy: compute every point order in the tiny ambient group and take their lcm"
        );
        assert!(
            describe_exhaustive_group_exponent_report(&exact_exponent)
                .contains("Exhaustive group exponent")
        );
        assert!(
            describe_exhaustive_group_exponent_report(&exact_exponent)
                .contains("exact exponent: 6")
        );
    }

    #[test]
    fn exponent_lower_bound_group_order_verification_visualization_stays_honest_about_scope() {
        let curve = crate::elliptic_curves::ShortWeierstrassCurve::<Fp<5>>::new(
            Fp::<5>::from_i64(0),
            Fp::<5>::from_i64(1),
        )
        .expect("valid curve");
        let sampled_point = curve
            .point(Fp::<5>::from_i64(2), Fp::<5>::from_i64(2))
            .expect("point should lie on the curve");
        let point_index = curve
            .points()
            .iter()
            .position(|candidate| candidate == &sampled_point)
            .expect("sample point should appear in the enumerated group");
        let mut sampler =
            move |upper_bound: usize| (point_index < upper_bound).then_some(point_index);

        let report = curve
            .group_exponent_by(
                crate::elliptic_curves::GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy: crate::elliptic_curves::PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("random-point exponent accumulation should succeed");
        let crate::elliptic_curves::GroupExponentReport::RandomPoints(accumulation) = report else {
            panic!("expected accumulation report");
        };
        let verification = curve
            .verify_exponent_lower_bound_by_group_order(
                &accumulation,
                crate::elliptic_curves::GroupOrderStrategy::Auto,
            )
            .expect("verification should succeed");

        assert_eq!(
            format_exponent_lower_bound_group_order_verification(&verification),
            "group order verifies #E(F_q) = 6 from lower bound 6"
        );

        let description = describe_exponent_lower_bound_group_order_verification(&verification);
        assert!(description.contains("Exponent lower-bound verification by group order"));
        assert!(description.contains("exponent lower bound: 6"));
        assert!(description.contains("verified group order: 6"));
        assert!(description.contains("does not by itself certify the exponent"));
    }

    #[test]
    fn visualizable_trait_is_hooked_up_for_curves_and_points() {
        let curve = f7_curve();
        let point = f7_point(2, 1);

        assert!(curve.describe().contains("Short-Weierstrass curve"));
        assert_eq!(point.format_compact(), format_point_compact(&point));
    }

    #[test]
    fn curve_display_works_over_q_too() {
        let curve = crate::elliptic_curves::ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1))
            .expect("valid curve");

        assert_eq!(curve.to_equation_string(), "y^2 = x^3 + (-1)x + (0)");
        assert_eq!(format!("{curve}"), curve.to_equation_string());
        assert_eq!(format_curve(&curve), "y^2 = x^3 + (-1)x");
    }
}
