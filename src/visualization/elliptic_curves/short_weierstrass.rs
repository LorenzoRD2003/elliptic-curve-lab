use core::fmt;

use num_bigint::BigUint;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    frobenius::group_order::GroupOrderRoute,
    short_weierstrass::{
        group_exponent::{
            ExponentAccumulationReport, ExponentAccumulationStep,
            ExponentLowerBoundGroupOrderVerification, GroupExponentReport, GroupExponentStrategy,
        },
        point_order::{
            ExhaustivePointOrderReport, HasseIntervalPointOrderReport,
            PointOrderFromMultipleReport, PointOrderReductionStep, PointOrderReport,
            PointOrderStrategyKind,
        },
        rational_torsion::{
            RationalTorsionGroupShape, RationalTorsionReport, RationalTorsionStrategy,
        },
    },
    traits::{
        CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure, FiniteGroupCurveModel,
        GroupCurveModel,
    },
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::frobenius::{
        describe_group_order_report, describe_hasse_multiple_search_report, format_hasse_interval,
        format_hasse_multiple_search_report,
    },
    shared::{format_field_elem as format_elem, parenthesize_if_needed, yes_no},
};

fn equation_string<F: Field>(curve: &ShortWeierstrassCurve<F>) -> String
where
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
pub(crate) fn format_curve<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    equation_string(curve)
}

/// Formats an affine point compactly.
pub(crate) fn format_point<F: Field>(point: &AffinePoint<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    point.to_coordinates_string()
}

/// Formats an affine point using the compact field-element visualization.
pub(crate) fn format_point_compact<F: Field>(point: &AffinePoint<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    match point {
        AffinePoint::Infinity => "O".to_string(),
        AffinePoint::Finite { x, y } => {
            format!("({}, {})", format_elem::<F>(x), format_elem::<F>(y))
        }
    }
}

impl<F: Field> Visualizable for AffinePoint<F>
where
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

impl<F: Field> Visualizable for ShortWeierstrassCurve<F>
where
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
fn describe_curve<F: Field>(curve: &ShortWeierstrassCurve<F>) -> String
where
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
fn describe_point<F: Field>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display,
{
    let mut lines = vec![
        "Curve point".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("point: {}", format_point(point)),
        format!("identity: {}", yes_no(curve.is_identity(point))),
        format!("on curve: {}", yes_no(curve.contains(point))),
    ];

    if let AffinePoint::Finite { x, y } = point {
        lines.push(format!("x-coordinate: {}", format_elem::<F>(x)));
        lines.push(format!("y-coordinate: {}", format_elem::<F>(y)));
    }

    lines.join("\n")
}

/// Formats a Mazur-shape rational-torsion classification over `Q`.
fn format_rational_torsion_group_shape(shape: RationalTorsionGroupShape) -> String {
    match shape {
        RationalTorsionGroupShape::Trivial => "{O}".to_string(),
        RationalTorsionGroupShape::Cyclic { order } => format!("ℤ/{order}ℤ"),
        RationalTorsionGroupShape::ProductZ2Z2m { m } => format!("ℤ/2ℤ × ℤ/{}ℤ", 2 * m),
    }
}

/// Describes the exact rational-torsion computation for a short-Weierstrass
/// curve over `Q`.
fn describe_rational_torsion_report(report: &RationalTorsionReport) -> String {
    let mut lines = vec![
        "Rational torsion over Q".to_string(),
        format!("source curve: {}", format_curve(report.original_curve())),
        format!("integral model: {}", format_curve(report.integral_model())),
        format!(
            "integral scale u: {}",
            format_elem::<crate::fields::Q>(report.scale())
        ),
    ];

    if report.original_curve() == report.integral_model() {
        lines.push("integral transport: source curve was already integral".to_string());
    } else {
        lines.push("integral transport: source curve was scaled before torsion search".to_string());
    }

    lines.push(format!(
        "strategy: {}",
        format_rational_torsion_strategy(report.strategy())
    ));
    lines.extend([
        format!(
            "group: {}",
            format_rational_torsion_group_shape(report.group().shape())
        ),
        format!("group cardinality: {}", report.group().cardinality()),
    ]);
    if let Some(candidate_count) = report.lutz_nagell_candidate_count() {
        lines.push(format!(
            "Lutz-Nagell candidates checked: {} ({} rejected)",
            candidate_count,
            report.lutz_nagell_rejected_candidate_count().unwrap_or(0)
        ));
    }
    lines.push("torsion points:".to_string());

    for point in report.points() {
        lines.push(format!("  {}", format_point_compact(point)));
    }

    lines.join("\n")
}

impl Visualizable for RationalTorsionReport {
    fn format_compact(&self) -> String {
        format!(
            "rational torsion: {}, {} point(s)",
            format_rational_torsion_group_shape(self.group().shape()),
            self.points().len()
        )
    }

    fn describe(&self) -> String {
        describe_rational_torsion_report(self)
    }
}

fn format_rational_torsion_strategy(strategy: RationalTorsionStrategy) -> &'static str {
    match strategy {
        RationalTorsionStrategy::LutzNagell => {
            "integral model -> Lutz-Nagell candidates -> Mazur-order verification"
        }
        RationalTorsionStrategy::GoodReductionHensel => {
            "integral model -> good reduction -> division-polynomial x-criteria -> Hensel -> Mazur-order verification"
        }
    }
}

/// Explains curve membership by comparing `y^2` with `x^3 + ax + b`.
fn describe_membership<F: Field>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
where
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
fn explain_add<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    left: &AffinePoint<F>,
    right: &AffinePoint<F>,
) -> Result<String, CurveError>
where
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
fn list_points<F>(curve: &ShortWeierstrassCurve<F>) -> String
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
fn describe_point_order<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
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

impl Visualizable for FiniteAbelianGroupStructure {
    fn format_compact(&self) -> String {
        format_invariant_factor_surface(*self)
    }

    fn describe(&self) -> String {
        [
            "Finite abelian group structure".to_string(),
            format!("group order: {}", self.order),
            format!("cyclic: {}", yes_no(self.cyclic)),
            format!("exponent: {}", self.exponent),
            format!(
                "invariant factors: {}",
                format_invariant_factor_surface(*self)
            ),
        ]
        .join("\n")
    }
}

/// Describes the finite abelian group structure of a small enumerated curve.
fn describe_group_structure<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let structure = curve.group_structure();

    [
        "Finite curve group structure".to_string(),
        format!("curve: {}", equation_string(curve)),
        format!("group order: {}", structure.order),
        format!("cyclic: {}", yes_no(structure.cyclic)),
        format!("exponent: {}", structure.exponent),
        format!(
            "invariant factors: {}",
            format_invariant_factor_surface(structure)
        ),
    ]
    .join("\n")
}

/// Returns a compact educational summary of the finite group structure.
fn summarize_group_structure<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: EnumerableFiniteField + SqrtField,
    F::Elem: VisualizableField + fmt::Display,
{
    let structure = curve.group_structure();

    [
        format!("cyclic: {}", yes_no(structure.cyclic)),
        format!("exponent: {}", structure.exponent),
        format!(
            "invariant factors: {}",
            format_invariant_factor_surface(structure)
        ),
    ]
    .join("\n")
}

/// Describes how many points have each exact order on a small finite curve.
fn describe_order_distribution<F>(curve: &ShortWeierstrassCurve<F>) -> String
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
fn summarize_order_distribution<F>(curve: &ShortWeierstrassCurve<F>) -> String
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
fn describe_scalar_mul<F>(
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
fn explain_point_order<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> String
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

mod group_exponent;
mod point_order;

use group_exponent::*;
use point_order::*;
pub(crate) use point_order::{
    describe_point_order_from_multiple_report, format_point_order_from_multiple_report,
};

#[cfg(test)]
use group_exponent::{
    describe_exhaustive_group_exponent_report,
    describe_exponent_lower_bound_group_order_verification, describe_group_exponent_report,
    format_exponent_lower_bound_group_order_verification, format_group_exponent_report,
};
#[cfg(test)]
use point_order::{describe_exhaustive_point_order_report, format_point_order_report};

#[cfg(test)]
mod tests;
