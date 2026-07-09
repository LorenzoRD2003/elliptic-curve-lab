use crate::isogenies::{
    class_group_action::CraterWalkReport,
    graphs::{IsogenyGraphNodeId, endomorphisms::CraterShape},
};
use crate::visualization::Visualizable;

/// Explains one deterministic crater walk labeled by a prime-norm ideal.
fn explain_crater_walk_report(report: &CraterWalkReport) -> String {
    let mut lines = vec![
        "Crater walk labeled by an ideal".to_string(),
        "--------------------------------".to_string(),
        format!("ideal norm ℓ: {}", report.ideal().norm()),
        format!("ideal root mod ℓ: {}", report.ideal().root_mod_ell()),
        format!(
            "crater shape: {}",
            format_crater_shape(report.crater_shape())
        ),
        format!("start: v{}", report.start().0),
        format!("visited: {}", format_node_path(report.visited())),
        format!(
            "cycle length: {}",
            report
                .cycle_length()
                .map(|length| length.to_string())
                .unwrap_or_else(|| "not closed".to_string())
        ),
    ];

    lines.push(
        "The walk follows a deterministic graph order on certified crater edges; it is not an oriented action of 𝔭."
            .to_string(),
    );

    lines.join("\n")
}

impl Visualizable for CraterWalkReport {
    fn format_compact(&self) -> String {
        match self.cycle_length() {
            Some(length) => format!(
                "crater walk: {} (cycle length {length})",
                format_node_path(self.visited())
            ),
            None => format!(
                "crater walk: {} (not closed)",
                format_node_path(self.visited())
            ),
        }
    }

    fn describe(&self) -> String {
        explain_crater_walk_report(self)
    }
}

fn format_crater_shape(shape: CraterShape) -> String {
    match shape {
        CraterShape::EmptyCertifiedCrater => "empty certified crater".to_string(),
        CraterShape::Singleton { self_loop_count } => {
            format!("singleton crater with {self_loop_count} self-loop edge(s)")
        }
        CraterShape::TwoVertex {
            directed_edge_count,
        } => {
            format!("two-vertex crater with {directed_edge_count} directed edge(s)")
        }
        CraterShape::Cycle { length } => format!("cycle of length {length}"),
        CraterShape::PartialOrAmbiguous => "partial or ambiguous crater".to_string(),
    }
}

fn format_node_path(nodes: &[IsogenyGraphNodeId]) -> String {
    if nodes.is_empty() {
        return "none".to_string();
    }

    nodes
        .iter()
        .map(|node| format!("v{}", node.0))
        .collect::<Vec<_>>()
        .join(" -> ")
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::explain_crater_walk_report;
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        endomorphisms::{
            quadratic_ideals::PrimeNormIdeal,
            quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
        },
    };
    use crate::fields::Fp7;
    use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};
    use crate::visualization::Visualizable;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    fn order_minus_23() -> ImaginaryQuadraticOrder {
        ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(1))
            .expect("D = -23 should define an imaginary quadratic maximal order")
    }

    fn f7_curve() -> ShortWeierstrassCurve<Fp7> {
        ShortWeierstrassCurve::<Fp7>::new(Fp7::from_i64(2), Fp7::from_i64(3))
            .expect("valid F_7 curve")
    }

    fn split_three_ideal() -> PrimeNormIdeal {
        PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
            .expect("3 splits in the order of discriminant -23")
    }

    #[test]
    fn crater_walk_explanation_mentions_path_and_context() {
        let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let report = graph
            .crater_walk_report(split_three_ideal(), IsogenyGraphNodeId(0))
            .expect("walk report should build");

        let explanation = explain_crater_walk_report(&report);

        assert!(explanation.contains("ideal norm ℓ: 3"));
        assert!(explanation.contains("ideal root mod ℓ: 1"));
        assert!(explanation.contains("two-vertex crater"));
        assert!(explanation.contains("visited: v0 -> v1 -> v0"));
        assert!(explanation.contains("cycle length: 2"));
        assert!(explanation.contains("deterministic graph order"));
        assert!(report.format_compact().contains("cycle length 2"));
    }
}
