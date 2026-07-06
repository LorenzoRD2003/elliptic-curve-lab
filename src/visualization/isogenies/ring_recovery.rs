use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, FrobeniusTraceCurveModel};
use crate::fields::traits::{EnumerableFiniteField, Field, FiniteField, SqrtField};
use crate::isogenies::graphs::endomorphisms::{
    EndomorphismRingLevelRecoveryError, EndomorphismRingLevelRecoveryReport,
    LocalEndomorphismRingLevelReport,
};
use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};
use crate::visualization::Visualizable;

/// Which global assembly to show in an educational ring-recovery walkthrough.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndomorphismRingRecoveryAssembly<'a> {
    /// Assemble a partial report from the first `count` primes in the supplied
    /// prime list.
    FirstPrimes {
        /// Label printed immediately before the assembled report.
        label: &'a str,
        /// Number of leading local reports to assemble.
        count: usize,
    },
    /// Assemble the complete report from every supplied local prime.
    Complete {
        /// Label printed immediately before the assembled report.
        label: &'a str,
    },
}

impl<'a> EndomorphismRingRecoveryAssembly<'a> {
    /// Requests a partial assembly from the first `count` local primes.
    pub fn first_primes(label: &'a str, count: usize) -> Self {
        Self::FirstPrimes { label, count }
    }

    /// Requests the complete assembly from all supplied local primes.
    pub fn complete(label: &'a str) -> Self {
        Self::Complete { label }
    }
}

/// Errors produced by the high-level educational ring-recovery walkthrough.
#[derive(Debug)]
pub enum EndomorphismRingRecoveryWalkthroughError {
    /// The walkthrough needs at least one local prime.
    EmptyPrimeList,
    /// A requested local section refers to a prime that was not recovered.
    MissingLocalPrime { prime: usize },
    /// A partial assembly requested more local reports than were recovered.
    AssemblyPrefixTooLong {
        count: usize,
        recovered_primes: usize,
    },
    /// The underlying graph-side recovery failed.
    Recovery(EndomorphismRingLevelRecoveryError),
}

impl fmt::Display for EndomorphismRingRecoveryWalkthroughError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyPrimeList => write!(
                formatter,
                "endomorphism-ring recovery walkthrough needs at least one local prime"
            ),
            Self::MissingLocalPrime { prime } => write!(
                formatter,
                "walkthrough requested a local section for ℓ = {prime}, but that prime was not recovered"
            ),
            Self::AssemblyPrefixTooLong {
                count,
                recovered_primes,
            } => write!(
                formatter,
                "walkthrough requested a partial assembly from {count} local reports, but only {recovered_primes} primes were recovered"
            ),
            Self::Recovery(error) => write!(formatter, "{error}"),
        }
    }
}

impl std::error::Error for EndomorphismRingRecoveryWalkthroughError {}

impl From<EndomorphismRingLevelRecoveryError> for EndomorphismRingRecoveryWalkthroughError {
    fn from(error: EndomorphismRingLevelRecoveryError) -> Self {
        Self::Recovery(error)
    }
}

/// Explains one local endomorphism-ring level recovery report.
///
/// The explanation emphasizes the Sutherland §3.3 identity
/// `v_ℓ(u) = e - δ`: the Frobenius conductor valuation `e`, the certified
/// distance to the floor `δ`, and the recovered local conductor exponent.
pub fn explain_local_endomorphism_ring_level_report(
    report: &LocalEndomorphismRingLevelReport,
) -> String {
    vec![
        "Local endomorphism-ring level recovery".to_string(),
        "----------------------------------------".to_string(),
        format!("node: v{}", report.node_id().0),
        format!("prime ℓ: {}", report.prime()),
        format!(
            "Frobenius conductor valuation e = v_ℓ(v): {}",
            report.frobenius_conductor_valuation()
        ),
        format!(
            "certified floor distance δ = dist(E, V_d): {}",
            report.distance_to_floor()
        ),
        format!(
            "recovered local exponent d = e - δ = v_ℓ(u): {}",
            report.recovered_conductor_valuation()
        ),
        format!("floor path: {}", format_floor_path(report)),
        "This recovers only one local component of End(E) ≅ O_u.".to_string(),
    ]
    .join("\n")
}

/// Explains a global endomorphism-ring level recovery report assembled from
/// local volcano evidence.
///
/// Complete reports identify the candidate order `O_u`; partial reports show
/// which prime factors of the Frobenius conductor `v` still lack local
/// recovery evidence.
pub fn explain_endomorphism_ring_level_recovery_report(
    report: &EndomorphismRingLevelRecoveryReport,
) -> String {
    let mut lines = vec![
        "Endomorphism-ring level recovery from volcano floors".to_string(),
        "-----------------------------------------------------".to_string(),
        format!(
            "node: {}",
            report
                .node_id()
                .map(|node| format!("v{}", node.0))
                .unwrap_or_else(|| "none (v = 1 needs no local reports)".to_string())
        ),
        format!(
            "Frobenius conductor v: {}",
            report.candidate_set().frobenius_conductor()
        ),
        format!(
            "fundamental discriminant D_K: {}",
            report.candidate_set().fundamental_discriminant().value()
        ),
        format!("local reports: {}", report.local_reports().len()),
        format!(
            "complete local coverage: {}",
            if report.is_complete() { "yes" } else { "no" }
        ),
        "The global conductor is assembled as u = ∏ℓ^{d_ℓ}.".to_string(),
        String::new(),
        "Local exponents:".to_string(),
    ];

    if report.local_reports().is_empty() {
        lines.push("  none".to_string());
    } else {
        lines.extend(
            report
                .local_reports()
                .iter()
                .map(format_local_exponent_line),
        );
    }

    lines.push(String::new());
    if report.is_complete() {
        lines.push(format!(
            "recovered conductor u: {}",
            report
                .recovered_conductor()
                .expect("complete reports expose a conductor")
        ));
        lines.push(format!(
            "recovered order: {}",
            format_order_label(
                report
                    .recovered_order()
                    .expect("complete reports expose an order")
            )
        ));
    } else {
        lines.push(format!(
            "missing primes: {}",
            report
                .missing_primes()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ));
        lines.push("recovered order: unavailable until every ℓ | v is covered".to_string());
    }

    lines.push(
        "This report assembles certified local data; it does not build or search additional graphs."
            .to_string(),
    );

    lines.join("\n")
}

/// Explains a small endomorphism-ring recovery walkthrough.
///
/// This helper is intentionally presentation-only: callers still build the
/// local and global recovery reports through the graph/endomorphism APIs.
/// The visualization layer only owns the educational layout that places
/// introductory text, local `ℓ`-volcano probes, and assembled global reports
/// in one readable block.
pub fn explain_endomorphism_ring_level_recovery_walkthrough(
    title: &str,
    introduction: &[&str],
    local_reports: &[(&str, &LocalEndomorphismRingLevelReport)],
    global_reports: &[(&str, &EndomorphismRingLevelRecoveryReport)],
) -> String {
    let mut lines = vec![title.to_string(), "-".repeat(title.chars().count())];
    lines.extend(introduction.iter().map(|line| (*line).to_string()));

    if !local_reports.is_empty() {
        lines.push(String::new());
        for (index, (label, report)) in local_reports.iter().enumerate() {
            if index > 0 {
                lines.push(String::new());
            }
            push_optional_label(&mut lines, label);
            lines.push(explain_local_endomorphism_ring_level_report(report));
        }
    }

    for (label, report) in global_reports {
        lines.push(String::new());
        push_optional_label(&mut lines, label);
        lines.push(explain_endomorphism_ring_level_recovery_report(report));
    }

    lines.join("\n")
}

/// Builds and explains a root-node ring-recovery walkthrough for a
/// short-Weierstrass curve.
///
/// This is a high-level educational convenience for runnable examples. It
/// builds the small `ℓ`-isogeny graphs for the supplied primes, recovers the
/// local root-node reports, assembles the requested partial or complete global
/// reports, and then delegates the text layout to
/// [`explain_endomorphism_ring_level_recovery_walkthrough`].
///
/// It is intentionally scoped to the current short-Weierstrass graph builder
/// rather than pretending to be a general graph orchestration framework.
pub fn explain_short_weierstrass_root_endomorphism_ring_level_recovery_walkthrough<F>(
    title: &str,
    introduction: &[&str],
    curve: ShortWeierstrassCurve<F>,
    primes: &[usize],
    local_sections: &[(&str, usize)],
    assemblies: &[EndomorphismRingRecoveryAssembly<'_>],
) -> Result<String, EndomorphismRingRecoveryWalkthroughError>
where
    F: Field + EnumerableFiniteField + SqrtField + FiniteField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
    <ShortWeierstrassCurve<F> as CurveModel>::BaseField: EnumerableFiniteField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + SqrtField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + FiniteField,
    <ShortWeierstrassCurve<F> as CurveModel>::Point: Clone + Eq + Hash + PartialEq,
    ShortWeierstrassCurve<F>: FrobeniusTraceCurveModel,
{
    let first_prime = primes
        .first()
        .copied()
        .ok_or(EndomorphismRingRecoveryWalkthroughError::EmptyPrimeList)?;
    let root = IsogenyGraphNodeId(0);
    let locals = primes
        .iter()
        .map(|&ell| recover_root_level_for_prime(curve.clone(), ell, root))
        .collect::<Result<Vec<_>, _>>()?;
    let candidate_set = candidate_set_for_root(curve, first_prime, root)?;

    let local_views = local_sections
        .iter()
        .map(|(label, prime)| {
            local_report_for_prime(&locals, *prime)
                .map(|report| (*label, report))
                .ok_or(
                    EndomorphismRingRecoveryWalkthroughError::MissingLocalPrime { prime: *prime },
                )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let globals = assemblies
        .iter()
        .map(|assembly| global_report_for_assembly(candidate_set.clone(), &locals, *assembly))
        .collect::<Result<Vec<_>, _>>()?;
    let global_views = assemblies
        .iter()
        .zip(globals.iter())
        .map(|(assembly, report)| (assembly.label(), report))
        .collect::<Vec<_>>();

    Ok(explain_endomorphism_ring_level_recovery_walkthrough(
        title,
        introduction,
        &local_views,
        &global_views,
    ))
}

impl Visualizable for LocalEndomorphismRingLevelReport {
    fn format_compact(&self) -> String {
        format!(
            "local End(E) level at ℓ = {} for v{}: d = {}",
            self.prime(),
            self.node_id().0,
            self.recovered_conductor_valuation()
        )
    }

    fn describe(&self) -> String {
        explain_local_endomorphism_ring_level_report(self)
    }
}

impl Visualizable for EndomorphismRingLevelRecoveryReport {
    fn format_compact(&self) -> String {
        if let Some(order) = self.recovered_order() {
            format!(
                "endomorphism-ring recovery: {} from {} local reports",
                format_order_label(order),
                self.local_reports().len()
            )
        } else {
            format!(
                "partial endomorphism-ring recovery: missing ℓ = {}",
                self.missing_primes()
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        }
    }

    fn describe(&self) -> String {
        explain_endomorphism_ring_level_recovery_report(self)
    }
}

fn push_optional_label(lines: &mut Vec<String>, label: &str) {
    if !label.is_empty() {
        lines.push(label.to_string());
    }
}

fn recover_root_level_for_prime<F>(
    curve: ShortWeierstrassCurve<F>,
    ell: usize,
    root: IsogenyGraphNodeId,
) -> Result<LocalEndomorphismRingLevelReport, EndomorphismRingLevelRecoveryError>
where
    F: Field + EnumerableFiniteField + SqrtField + FiniteField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
    <ShortWeierstrassCurve<F> as CurveModel>::BaseField: EnumerableFiniteField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + SqrtField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + FiniteField,
    <ShortWeierstrassCurve<F> as CurveModel>::Point: Clone + Eq + Hash + PartialEq,
    ShortWeierstrassCurve<F>: FrobeniusTraceCurveModel,
{
    let graph = IsogenyGraphBuilder::new(curve, ell)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    Ok(graph.recover_endomorphism_ring_level_at(root, &BigUint::from(ell))?)
}

fn candidate_set_for_root<F>(
    curve: ShortWeierstrassCurve<F>,
    ell: usize,
    root: IsogenyGraphNodeId,
) -> Result<
    crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet,
    EndomorphismRingLevelRecoveryError,
>
where
    F: Field + EnumerableFiniteField + SqrtField + FiniteField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
    <ShortWeierstrassCurve<F> as CurveModel>::BaseField: EnumerableFiniteField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + SqrtField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + FiniteField,
    <ShortWeierstrassCurve<F> as CurveModel>::Point: Clone + Eq + Hash + PartialEq,
    ShortWeierstrassCurve<F>: FrobeniusTraceCurveModel,
{
    let graph = IsogenyGraphBuilder::new(curve, ell)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    Ok(graph.node_endomorphism_candidates(root)?)
}

fn local_report_for_prime(
    reports: &[LocalEndomorphismRingLevelReport],
    prime: usize,
) -> Option<&LocalEndomorphismRingLevelReport> {
    let prime = BigUint::from(prime);
    reports.iter().find(|report| report.prime() == &prime)
}

fn global_report_for_assembly(
    candidate_set: crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet,
    locals: &[LocalEndomorphismRingLevelReport],
    assembly: EndomorphismRingRecoveryAssembly<'_>,
) -> Result<EndomorphismRingLevelRecoveryReport, EndomorphismRingRecoveryWalkthroughError> {
    let reports = match assembly {
        EndomorphismRingRecoveryAssembly::FirstPrimes { count, .. } => {
            if count > locals.len() {
                return Err(
                    EndomorphismRingRecoveryWalkthroughError::AssemblyPrefixTooLong {
                        count,
                        recovered_primes: locals.len(),
                    },
                );
            }
            locals[..count].to_vec()
        }
        EndomorphismRingRecoveryAssembly::Complete { .. } => locals.to_vec(),
    };

    Ok(EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set,
        reports,
    )?)
}

impl EndomorphismRingRecoveryAssembly<'_> {
    fn label(&self) -> &str {
        match self {
            Self::FirstPrimes { label, .. } | Self::Complete { label } => label,
        }
    }
}

fn format_floor_path(report: &LocalEndomorphismRingLevelReport) -> String {
    report
        .floor_path()
        .path()
        .iter()
        .map(|node| format!("v{}", node.0))
        .collect::<Vec<_>>()
        .join(" -> ")
}

fn format_local_exponent_line(report: &LocalEndomorphismRingLevelReport) -> String {
    format!(
        "  ℓ = {}: e = {}, δ = {}, d = {}",
        report.prime(),
        report.frobenius_conductor_valuation(),
        report.distance_to_floor(),
        report.recovered_conductor_valuation()
    )
}

fn format_order_label(
    order: &crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder,
) -> String {
    format!(
        "O_{} (Δ = {})",
        order.conductor(),
        order.discriminant().value()
    )
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
    };
    use crate::fields::Fp17;
    use crate::isogenies::graphs::{
        IsogenyGraphNodeId,
        endomorphisms::{
            EndomorphismRingLevelRecoveryReport, LocalEndomorphismRingLevelReport,
            ShortestFloorPathReport,
        },
    };
    use crate::visualization::isogenies::{
        EndomorphismRingRecoveryAssembly, EndomorphismRingRecoveryWalkthroughError,
        explain_endomorphism_ring_level_recovery_report,
        explain_endomorphism_ring_level_recovery_walkthrough,
        explain_local_endomorphism_ring_level_report,
        explain_short_weierstrass_root_endomorphism_ring_level_recovery_walkthrough,
    };
    use crate::visualization::traits::Visualizable;

    fn candidate_set(discriminant: i64) -> EndomorphismRingCandidateSet {
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(discriminant))
            .expect("test discriminant should produce candidate orders")
    }

    fn local_report(
        prime: u8,
        discriminant: i64,
        path: Vec<IsogenyGraphNodeId>,
    ) -> LocalEndomorphismRingLevelReport {
        let prime = BigUint::from(prime);
        let local_view = candidate_set(discriminant)
            .local_view_at(&prime)
            .expect("test prime should be valid");
        let start = path.first().copied().expect("path should be nonempty");
        let floor = path.last().copied().expect("path should be nonempty");
        let floor_path = ShortestFloorPathReport::new(prime, start, floor, path);

        LocalEndomorphismRingLevelReport::from_local_view_and_floor_path(local_view, floor_path)
            .expect("test local report should be compatible")
    }

    #[test]
    fn local_recovery_explanation_mentions_e_delta_and_recovered_d() {
        let report = local_report(2, -64, vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]);

        let explanation = explain_local_endomorphism_ring_level_report(&report);

        assert!(explanation.contains("Local endomorphism-ring level recovery"));
        assert!(explanation.contains("prime ℓ: 2"));
        assert!(explanation.contains("e = v_ℓ(v): 2"));
        assert!(explanation.contains("δ = dist(E, V_d): 1"));
        assert!(explanation.contains("d = e - δ = v_ℓ(u): 1"));
        assert!(explanation.contains("floor path: v0 -> v1"));
        assert_eq!(
            report.format_compact(),
            "local End(E) level at ℓ = 2 for v0: d = 1"
        );
    }

    #[test]
    fn complete_global_recovery_explanation_mentions_recovered_order() {
        let report = EndomorphismRingLevelRecoveryReport::from_local_reports(
            candidate_set(-108),
            vec![
                local_report(2, -108, vec![IsogenyGraphNodeId(0)]),
                local_report(3, -108, vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(2)]),
            ],
        )
        .expect("complete report should assemble");

        let explanation = explain_endomorphism_ring_level_recovery_report(&report);

        assert!(explanation.contains("Endomorphism-ring level recovery from volcano floors"));
        assert!(explanation.contains("Frobenius conductor v: 6"));
        assert!(explanation.contains("complete local coverage: yes"));
        assert!(explanation.contains("u = ∏ℓ^{d_ℓ}"));
        assert!(explanation.contains("ℓ = 2: e = 1, δ = 0, d = 1"));
        assert!(explanation.contains("ℓ = 3: e = 1, δ = 1, d = 0"));
        assert!(explanation.contains("recovered conductor u: 2"));
        assert!(explanation.contains("recovered order: O_2"));
        assert!(report.format_compact().contains("O_2"));
    }

    #[test]
    fn partial_global_recovery_explanation_mentions_missing_primes() {
        let report = EndomorphismRingLevelRecoveryReport::from_local_reports(
            candidate_set(-108),
            vec![local_report(2, -108, vec![IsogenyGraphNodeId(0)])],
        )
        .expect("partial report should assemble");

        let explanation = explain_endomorphism_ring_level_recovery_report(&report);

        assert!(explanation.contains("complete local coverage: no"));
        assert!(explanation.contains("missing primes: 3"));
        assert!(explanation.contains("recovered order: unavailable"));
        assert!(report.format_compact().contains("missing ℓ = 3"));
    }

    #[test]
    fn recovery_walkthrough_places_intro_local_and_global_reports() {
        let local = local_report(3, -108, vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(2)]);
        let global = EndomorphismRingLevelRecoveryReport::from_local_reports(
            candidate_set(-108),
            vec![local.clone()],
        )
        .expect("partial report should assemble");

        let explanation = explain_endomorphism_ring_level_recovery_walkthrough(
            "Tiny walkthrough",
            &["intro line"],
            &[("Local probe:", &local)],
            &[("Partial assembly:", &global)],
        );

        assert!(explanation.contains("Tiny walkthrough\n----------------"));
        assert!(explanation.contains("intro line"));
        assert!(explanation.contains("Local probe:"));
        assert!(explanation.contains("prime ℓ: 3"));
        assert!(explanation.contains("Partial assembly:"));
        assert!(explanation.contains("missing primes: 2"));
    }

    #[test]
    fn short_weierstrass_root_walkthrough_builds_reports_for_examples() {
        let curve =
            ShortWeierstrassCurve::<Fp17>::new(Fp17::from_i64(11), Fp17::from_i64(5)).unwrap();

        let explanation =
            explain_short_weierstrass_root_endomorphism_ring_level_recovery_walkthrough(
                "Root walkthrough",
                &["computed from a curve"],
                curve,
                &[2],
                &[("Local 2-volcano probe:", 2)],
                &[EndomorphismRingRecoveryAssembly::complete(
                    "Complete assembly:",
                )],
            )
            .expect("walkthrough should build from the tiny graph");

        assert!(explanation.contains("Root walkthrough"));
        assert!(explanation.contains("computed from a curve"));
        assert!(explanation.contains("Local 2-volcano probe:"));
        assert!(explanation.contains("prime ℓ: 2"));
        assert!(explanation.contains("Complete assembly:"));
        assert!(explanation.contains("recovered order:"));
    }

    #[test]
    fn short_weierstrass_root_walkthrough_rejects_missing_local_prime_sections() {
        let curve =
            ShortWeierstrassCurve::<Fp17>::new(Fp17::from_i64(11), Fp17::from_i64(5)).unwrap();

        let error = explain_short_weierstrass_root_endomorphism_ring_level_recovery_walkthrough(
            "Root walkthrough",
            &[],
            curve,
            &[2],
            &[("Missing 3-volcano probe:", 3)],
            &[],
        )
        .expect_err("local sections must refer to recovered primes");

        assert!(matches!(
            error,
            EndomorphismRingRecoveryWalkthroughError::MissingLocalPrime { prime: 3 }
        ));
    }
}
