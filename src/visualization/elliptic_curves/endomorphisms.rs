use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet, quadratic_ideals::QuadraticPrimeBehavior,
    quadratic_orders::QuadraticOrderCoverRelation,
};
use crate::visualization::Visualizable;
use crate::visualization::shared::{comma_list, format_order_conductor_label};

fn format_cover_relation(relation: &QuadraticOrderCoverRelation) -> String {
    format!(
        "{} -> {} [index {}]",
        format_order_conductor_label(relation.overorder().conductor()),
        format_order_conductor_label(relation.suborder().conductor()),
        relation.index()
    )
}

/// Describes the Hasse-diagram view of the Frobenius-compatible candidate orders.
fn describe_endomorphism_ring_candidate_poset(
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
            format_order_conductor_label(candidate_set.frobenius_order().conductor())
        ),
        format!(
            "maximal order: {} = O_K",
            format_order_conductor_label(candidate_set.maximal_order().conductor())
        ),
        format!(
            "candidate orders: {}",
            comma_list(
                candidate_set
                    .candidate_orders()
                    .iter()
                    .map(|order| format_order_conductor_label(order.conductor())),
            )
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

fn describe_quadratic_prime_behavior(behavior: &QuadraticPrimeBehavior) -> String {
    let mut lines = vec![
        "Quadratic prime behavior".to_string(),
        "------------------------".to_string(),
        format!("status: {}", format_quadratic_prime_behavior(behavior)),
    ];

    match behavior {
        QuadraticPrimeBehavior::Split { roots } => {
            lines.push(format!("roots mod ℓ: {}, {}", roots.0, roots.1));
            lines.push("horizontal interpretation: possible crater directions".to_string());
        }
        QuadraticPrimeBehavior::Inert => {
            lines.push("horizontal interpretation: no horizontal ℓ-isogeny direction".to_string());
        }
        QuadraticPrimeBehavior::Ramified { root } => {
            lines.push(format!("repeated root mod ℓ: {root}"));
            lines.push("horizontal interpretation: degenerate local direction".to_string());
        }
        QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor => {
            lines.push("horizontal interpretation: ℓ is not invertible in this order".to_string());
        }
    }

    lines.join("\n")
}

fn format_quadratic_prime_behavior(behavior: &QuadraticPrimeBehavior) -> &'static str {
    match behavior {
        QuadraticPrimeBehavior::Split { .. } => "split",
        QuadraticPrimeBehavior::Inert => "inert",
        QuadraticPrimeBehavior::Ramified { .. } => "ramified",
        QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor => {
            "non-invertible conductor prime"
        }
    }
}

impl Visualizable for EndomorphismRingCandidateSet {
    fn format_compact(&self) -> String {
        format!("candidate orders for Δ_π = {}", self.discriminant().value())
    }

    fn describe(&self) -> String {
        describe_endomorphism_ring_candidate_poset(self)
    }
}

impl Visualizable for QuadraticPrimeBehavior {
    fn format_compact(&self) -> String {
        format!(
            "quadratic prime behavior: {}",
            format_quadratic_prime_behavior(self)
        )
    }

    fn describe(&self) -> String {
        describe_quadratic_prime_behavior(self)
    }
}

#[cfg(test)]
mod tests {
    use super::{describe_endomorphism_ring_candidate_poset, describe_quadratic_prime_behavior};
    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_ideals::QuadraticPrimeBehavior,
        quadratic_orders::QuadraticDiscriminant,
    };
    use crate::visualization::Visualizable;

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

    #[test]
    fn quadratic_prime_behavior_description_reports_horizontal_interpretation() {
        let split = QuadraticPrimeBehavior::Split {
            roots: (1u8.into(), 2u8.into()),
        };
        let inert = QuadraticPrimeBehavior::Inert;
        let ramified = QuadraticPrimeBehavior::Ramified { root: 0u8.into() };
        let non_invertible = QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor;

        assert!(describe_quadratic_prime_behavior(&split).contains("possible crater directions"));
        assert!(split.format_compact().contains("split"));
        assert!(
            inert
                .describe()
                .contains("no horizontal ℓ-isogeny direction")
        );
        assert!(ramified.describe().contains("degenerate local direction"));
        assert!(non_invertible.describe().contains("not invertible"));
    }
}
