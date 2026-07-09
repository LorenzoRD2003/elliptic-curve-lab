use crate::isogenies::graphs::endomorphisms::{
    EndomorphismRingLevelRecoveryReport, LocalEndomorphismRingLevelReport,
};
use crate::visualization::{
    Visualizable,
    shared::{comma_list, format_imaginary_quadratic_order_label as format_order_label, yes_no},
};

/// Explains one local endomorphism-ring level recovery report.
///
/// The explanation emphasizes the Sutherland §3.3 identity
/// `v_ℓ(u) = e - δ`: the Frobenius conductor valuation `e`, the certified
/// distance to the floor `δ`, and the recovered local conductor exponent.
fn explain_local_endomorphism_ring_level_report(
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
fn explain_endomorphism_ring_level_recovery_report(
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
        format!("complete local coverage: {}", yes_no(report.is_complete())),
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
            comma_list(report.missing_primes().iter().map(ToString::to_string))
        ));
        lines.push("recovered order: unavailable until every ℓ | v is covered".to_string());
    }

    lines.push(
        "This report assembles certified local data produced by the graph recovery route."
            .to_string(),
    );

    lines.join("\n")
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
                comma_list(self.missing_primes().iter().map(ToString::to_string))
            )
        }
    }

    fn describe(&self) -> String {
        explain_endomorphism_ring_level_recovery_report(self)
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

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::*;
    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
    };
    use crate::isogenies::graphs::{
        IsogenyGraphNodeId,
        endomorphisms::{
            EndomorphismRingLevelRecoveryReport, LocalEndomorphismRingLevelReport,
            ShortestFloorPathReport,
        },
    };
    use crate::visualization::Visualizable;

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
}
