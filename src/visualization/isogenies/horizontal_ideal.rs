use crate::isogenies::class_group_action::{
    HorizontalIdealReport, HorizontalIdealStatus, HorizontalIdealWitness,
};
use crate::visualization::Visualizable;

/// Explains one horizontal ideal compatibility report.
///
/// The explanation stays deliberately modest: it describes compatibility
/// between existing crater evidence and a caller-supplied prime-norm ideal, not
/// a computed class-group action.
fn explain_horizontal_ideal_report(report: &HorizontalIdealReport) -> String {
    let mut lines = vec![
        "Horizontal ideal compatibility".to_string(),
        "------------------------------".to_string(),
        format!(
            "status: {}",
            format_horizontal_ideal_status(report.status())
        ),
    ];

    match report.witness() {
        Some(witness) => lines.extend(format_witness_lines(witness)),
        None => lines.push("witness: none".to_string()),
    }

    lines.push(
        "This is compatibility evidence only; no ideal action or orientation is claimed."
            .to_string(),
    );

    lines.join("\n")
}

/// Explains a list of horizontal ideal compatibility reports.
fn explain_horizontal_ideal_reports(reports: &[HorizontalIdealReport]) -> String {
    let compatible = reports
        .iter()
        .filter(|report| report.status() == HorizontalIdealStatus::CertifiedCompatible)
        .count();
    let uncertified = reports
        .iter()
        .filter(|report| report.status() == HorizontalIdealStatus::EdgeNotCertifiedHorizontal)
        .count();
    let mismatched = reports
        .iter()
        .filter(|report| report.status() == HorizontalIdealStatus::DegreeMismatch)
        .count();

    let mut lines = vec![
        "Horizontal ideal compatibility reports".to_string(),
        "--------------------------------------".to_string(),
        format!("reports: {}", reports.len()),
        format!("certified compatible: {compatible}"),
        format!("edge not certified horizontal: {uncertified}"),
        format!("degree/norm mismatch: {mismatched}"),
        String::new(),
        "Entries:".to_string(),
    ];

    if reports.is_empty() {
        lines.push("  none".to_string());
    } else {
        lines.extend(reports.iter().map(format_report_line));
    }

    lines.push(String::new());
    lines.push(
        "These reports annotate crater evidence with supplied ideals; they do not compute [𝔞] * E."
            .to_string(),
    );

    lines.join("\n")
}

impl Visualizable for HorizontalIdealReport {
    fn format_compact(&self) -> String {
        match self.witness() {
            Some(witness) => format!(
                "horizontal ideal witness: e{} at ℓ = {}",
                witness.edge().edge_id().0,
                witness.prime()
            ),
            None => format!(
                "horizontal ideal report: {}",
                format_horizontal_ideal_status(self.status())
            ),
        }
    }

    fn describe(&self) -> String {
        explain_horizontal_ideal_report(self)
    }
}

fn format_witness_lines(witness: &HorizontalIdealWitness) -> Vec<String> {
    vec![
        format!("edge: e{}", witness.edge().edge_id().0),
        format!(
            "endpoints: v{} -> v{}",
            witness.edge().source().0,
            witness.edge().target().0
        ),
        format!("volcano prime ℓ: {}", witness.prime()),
        format!("ideal norm: {}", witness.ideal().norm()),
        format!("ideal root mod ℓ: {}", witness.ideal().root_mod_ell()),
    ]
}

fn format_report_line(report: &HorizontalIdealReport) -> String {
    match report.witness() {
        Some(witness) => format!(
            "  e{}: {} at ℓ = {}",
            witness.edge().edge_id().0,
            format_horizontal_ideal_status(report.status()),
            witness.prime()
        ),
        None => format!("  {}", format_horizontal_ideal_status(report.status())),
    }
}

fn format_horizontal_ideal_status(status: HorizontalIdealStatus) -> &'static str {
    match status {
        HorizontalIdealStatus::CertifiedCompatible => "certified compatible",
        HorizontalIdealStatus::EdgeNotCertifiedHorizontal => "edge not certified horizontal",
        HorizontalIdealStatus::DegreeMismatch => "degree/norm mismatch",
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::*;
    use crate::elliptic_curves::endomorphisms::{
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    };
    use crate::isogenies::{
        class_group_action::{HorizontalIdealReport, HorizontalIdealStatus},
        graphs::{
            IsogenyGraphEdgeId, IsogenyGraphNodeId,
            endomorphisms::{HorizontalEdgeReport, HorizontalEdgeStatus},
        },
    };
    use crate::visualization::Visualizable;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    fn order_minus_23() -> ImaginaryQuadraticOrder {
        ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(1))
            .expect("D = -23 should define an imaginary quadratic maximal order")
    }

    fn edge(status: HorizontalEdgeStatus) -> HorizontalEdgeReport {
        HorizontalEdgeReport::new(
            IsogenyGraphEdgeId(4),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2),
            status,
        )
    }

    fn ideal() -> PrimeNormIdeal {
        PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
            .expect("3 splits in the order of discriminant -23")
    }

    #[test]
    fn horizontal_ideal_report_explanation_mentions_witness_data() {
        let report = HorizontalIdealReport::from_certified_edge_and_ideal(
            edge(HorizontalEdgeStatus::CertifiedByAltitude),
            bu(3),
            ideal(),
        );

        let explanation = explain_horizontal_ideal_report(&report);

        assert!(explanation.contains("certified compatible"));
        assert!(explanation.contains("e4"));
        assert!(explanation.contains("volcano prime ℓ: 3"));
        assert!(explanation.contains("no ideal action"));
        assert!(report.format_compact().contains("e4"));
    }

    #[test]
    fn horizontal_ideal_reports_summary_counts_statuses() {
        let reports = vec![
            HorizontalIdealReport::from_certified_edge_and_ideal(
                edge(HorizontalEdgeStatus::CertifiedByAltitude),
                bu(3),
                ideal(),
            ),
            HorizontalIdealReport::from_certified_edge_and_ideal(
                edge(HorizontalEdgeStatus::SuspectedByWeakSurfaceEvidence),
                bu(3),
                ideal(),
            ),
            HorizontalIdealReport::from_certified_edge_and_ideal(
                edge(HorizontalEdgeStatus::CertifiedByAltitude),
                bu(5),
                ideal(),
            ),
        ];

        let explanation = explain_horizontal_ideal_reports(&reports);

        assert!(explanation.contains("reports: 3"));
        assert!(explanation.contains("certified compatible: 1"));
        assert!(explanation.contains("edge not certified horizontal: 1"));
        assert!(explanation.contains("degree/norm mismatch: 1"));
    }

    #[test]
    fn compact_report_mentions_non_witness_status() {
        let report = HorizontalIdealReport::from_certified_edge_and_ideal(
            edge(HorizontalEdgeStatus::CertifiedByAltitude),
            bu(5),
            ideal(),
        );

        assert_eq!(report.status(), HorizontalIdealStatus::DegreeMismatch);
        assert!(report.format_compact().contains("degree/norm mismatch"));
    }
}
