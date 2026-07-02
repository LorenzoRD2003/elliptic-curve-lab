use crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet;
use crate::elliptic_curves::endomorphisms::quadratic_orders::QuadraticOrderCoverRelation;
use crate::visualization::traits::Visualizable;

fn format_order_label(conductor: &num_bigint::BigUint) -> String {
    format!("O_{}", conductor)
}

fn format_cover_relation(relation: &QuadraticOrderCoverRelation) -> String {
    format!(
        "{} -> {} [index {}]",
        format_order_label(relation.overorder().conductor()),
        format_order_label(relation.suborder().conductor()),
        relation.index()
    )
}

/// Describes the Hasse-diagram view of the Frobenius-compatible candidate orders.
pub fn describe_endomorphism_ring_candidate_poset(
    candidate_set: &EndomorphismRingCandidateSet,
) -> String {
    let mut lines = vec![
        "Endomorphism-ring candidate poset".to_string(),
        format!(
            "Frobenius discriminant Δ_π: {}",
            candidate_set.discriminant().value()
        ),
        format!(
            "fundamental discriminant D_K: {}",
            candidate_set.fundamental_discriminant().value()
        ),
        format!(
            "Frobenius conductor v: {}",
            candidate_set.frobenius_conductor()
        ),
        format!(
            "Frobenius order: {} = ℤ[π]",
            format_order_label(candidate_set.frobenius_order().conductor())
        ),
        format!(
            "maximal order: {} = O_K",
            format_order_label(candidate_set.maximal_order().conductor())
        ),
        format!(
            "candidate orders: {}",
            candidate_set
                .candidate_orders()
                .iter()
                .map(|order| format_order_label(order.conductor()))
                .collect::<Vec<_>>()
                .join(", ")
        ),
        "Hasse cover relations (larger order -> immediately contained smaller order):".to_string(),
    ];

    let cover_relations = candidate_set.hasse_cover_relations();
    if cover_relations.is_empty() {
        lines.push("  none".to_string());
    } else {
        for relation in &cover_relations {
            lines.push(format!("  {}", format_cover_relation(relation)));
        }
    }

    lines.push(
        "interpretation: an edge O_g -> O_f means O_f ⊆ O_g and the label is the relative index [O_g : O_f]".to_string(),
    );

    lines.join("\n")
}

impl Visualizable for EndomorphismRingCandidateSet {
    fn format_compact(&self) -> String {
        format!("candidate orders for Δ_π = {}", self.discriminant().value())
    }

    fn describe(&self) -> String {
        describe_endomorphism_ring_candidate_poset(self)
    }
}

#[cfg(test)]
mod tests {

    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
    };
    use crate::visualization::elliptic_curves::describe_endomorphism_ring_candidate_poset;
    use crate::visualization::traits::Visualizable;

    #[test]
    fn endomorphism_candidate_poset_description_reports_hasse_edges_for_a_branching_example() {
        let candidate_set =
            EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-144))
                .expect("-144 should define a candidate-order poset with branching");

        let text = describe_endomorphism_ring_candidate_poset(&candidate_set);

        assert!(text.contains("Endomorphism-ring candidate poset"));
        assert!(text.contains("Frobenius conductor v: 6"));
        assert!(text.contains("candidate orders: O_1, O_2, O_3, O_6"));
        assert!(text.contains("O_1 -> O_2 [index 2]"));
        assert!(text.contains("O_1 -> O_3 [index 3]"));
        assert!(text.contains("O_2 -> O_6 [index 3]"));
        assert!(text.contains("O_3 -> O_6 [index 2]"));
        assert!(text.contains("the label is the relative index [O_g : O_f]"));
    }

    #[test]
    fn endomorphism_candidate_poset_is_visualizable() {
        let candidate_set =
            EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-48))
                .expect("-48 should define a candidate-order poset");

        assert_eq!(
            candidate_set.format_compact(),
            "candidate orders for Δ_π = -48"
        );
        assert!(candidate_set.describe().contains("O_1 -> O_2 [index 2]"));
    }
}
