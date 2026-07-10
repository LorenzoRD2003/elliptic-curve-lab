use num_bigint::BigInt;
use num_traits::Signed;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{
        BinaryQuadraticForm, QuadraticClassGroupCayleyTable, QuadraticClassGroupGeneratedSubgroup,
    },
    candidate_sets::EndomorphismRingCandidateSet,
    quadratic_ideals::QuadraticPrimeBehavior,
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

fn compact_binary_quadratic_form(form: &BinaryQuadraticForm) -> String {
    let (a, b, c) = form.coefficients();
    format!("({a},{b},{c})")
}

fn describe_binary_quadratic_form(form: &BinaryQuadraticForm) -> String {
    let mut lines = vec![
        "Binary quadratic form".to_string(),
        "---------------------".to_string(),
        format!("form: {}", compact_binary_quadratic_form(form)),
        format!("polynomial: {}", format_quadratic_polynomial(form)),
        format!("discriminant Δ = b² − 4ac: {}", form.discriminant()),
        format!("primitive: {}", yes_no(form.is_primitive())),
        format!("positive definite: {}", yes_no(form.is_positive_definite())),
        format!(
            "reduced positive definite: {}",
            yes_no(form.is_reduced_positive_definite())
        ),
    ];

    lines.push("coefficient storage: integral ternary (a,b,c)".to_string());
    lines.join("\n")
}

fn format_quadratic_polynomial(form: &BinaryQuadraticForm) -> String {
    let mut text = format!("{}x²", form.a());
    push_signed_term(&mut text, form.b(), "xy");
    push_signed_term(&mut text, form.c(), "y²");
    text
}

fn push_signed_term(text: &mut String, coefficient: &BigInt, monomial: &str) {
    if coefficient.is_negative() {
        text.push_str(&format!(" - {}{monomial}", -coefficient));
    } else {
        text.push_str(&format!(" + {coefficient}{monomial}"));
    }
}

fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

fn cayley_labels(table: &QuadraticClassGroupCayleyTable) -> Vec<String> {
    if table.class_number() == 0 {
        return Vec::new();
    }

    if let Some(labels) = cyclic_labels(table) {
        return labels;
    }

    if table.class_number() == 4 && all_nonidentity_elements_square_to_identity(table) {
        return vec![
            "e".to_string(),
            "a".to_string(),
            "b".to_string(),
            "ab".to_string(),
        ];
    }

    (0..table.class_number())
        .map(|index| {
            if index == 0 {
                "e".to_string()
            } else {
                format!("f{index}")
            }
        })
        .collect()
}

fn cyclic_labels(table: &QuadraticClassGroupCayleyTable) -> Option<Vec<String>> {
    for generator in 1..table.class_number() {
        let mut labels = vec![None; table.class_number()];
        labels[0] = Some("e".to_string());

        let mut current = 0usize;
        for exponent in 1..table.class_number() {
            current = table.product_index(current, generator)?;
            if labels[current].is_some() {
                break;
            }
            labels[current] = Some(power_label(exponent));
        }

        if table.product_index(current, generator)? == 0 && labels.iter().all(Option::is_some) {
            return labels.into_iter().collect();
        }
    }

    None
}

fn power_label(exponent: usize) -> String {
    if exponent == 1 {
        "g".to_string()
    } else {
        format!("g{}", superscript_number(exponent))
    }
}

fn superscript_number(value: usize) -> String {
    value
        .to_string()
        .chars()
        .map(|digit| match digit {
            '0' => '⁰',
            '1' => '¹',
            '2' => '²',
            '3' => '³',
            '4' => '⁴',
            '5' => '⁵',
            '6' => '⁶',
            '7' => '⁷',
            '8' => '⁸',
            '9' => '⁹',
            _ => digit,
        })
        .collect()
}

fn all_nonidentity_elements_square_to_identity(table: &QuadraticClassGroupCayleyTable) -> bool {
    (1..table.class_number()).all(|index| table.product_index(index, index) == Some(0))
}

fn describe_cayley_table(table: &QuadraticClassGroupCayleyTable) -> String {
    let labels = cayley_labels(table);
    let width = labels
        .iter()
        .map(|label| label.chars().count())
        .chain(["*".len()].into_iter())
        .max()
        .unwrap_or(1);
    let cell = |label: &str| format!("{label:>width$}");

    let mut lines = vec![
        "Quadratic class-group Cayley table".to_string(),
        "----------------------------------".to_string(),
        format!("discriminant D: {}", table.discriminant().value()),
        format!("class number h(D): {}", table.class_number()),
        "representatives:".to_string(),
    ];

    for (index, representative) in table.representatives().iter().enumerate() {
        lines.push(format!(
            "  {} = {}",
            labels[index],
            representative.format_compact()
        ));
    }

    lines.push("products:".to_string());
    let header = labels
        .iter()
        .map(|label| cell(label))
        .collect::<Vec<_>>()
        .join(" ");
    lines.push(format!("  {} | {}", cell("*"), header));
    lines.push(format!(
        "  {}-+-{}",
        "-".repeat(width),
        "-".repeat(header.len())
    ));

    for (row_index, row) in table.products().iter().enumerate() {
        let entries = row
            .iter()
            .map(|&product_index| cell(&labels[product_index]))
            .collect::<Vec<_>>()
            .join(" ");
        lines.push(format!("  {} | {}", cell(&labels[row_index]), entries));
    }

    lines.push("construction cost: Θ(h(D)²) class compositions".to_string());
    lines.join("\n")
}

fn describe_generated_subgroup(subgroup: &QuadraticClassGroupGeneratedSubgroup) -> String {
    let mut lines = vec![
        "Quadratic class-group generated subgroup".to_string(),
        "----------------------------------------".to_string(),
        format!("discriminant D: {}", subgroup.discriminant().value()),
        format!("generator g: {}", subgroup.generator().format_compact()),
        format!("subgroup order: {}", subgroup.order()),
        format!("ambient class number h(D): {}", subgroup.class_number()),
        format!(
            "generates whole class group: {}",
            yes_no(subgroup.is_whole_class_group())
        ),
        "powers:".to_string(),
    ];

    for (exponent, element) in subgroup.elements().iter().enumerate() {
        lines.push(format!(
            "  {} = {}",
            subgroup_power_label(exponent),
            element.format_compact()
        ));
    }

    lines.push(
        "interpretation: this is an algebraic cyclic subgroup of the form class group, not a certified crater action."
            .to_string(),
    );
    lines.push("construction cost: up to h(D) class compositions, O(h(D) · C) total".to_string());
    lines.join("\n")
}

fn subgroup_power_label(exponent: usize) -> String {
    if exponent == 0 {
        "g⁰".to_string()
    } else if exponent == 1 {
        "g¹".to_string()
    } else {
        format!("g{}", superscript_number(exponent))
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

impl Visualizable for BinaryQuadraticForm {
    fn format_compact(&self) -> String {
        compact_binary_quadratic_form(self)
    }

    fn describe(&self) -> String {
        describe_binary_quadratic_form(self)
    }
}

impl Visualizable for QuadraticClassGroupCayleyTable {
    fn format_compact(&self) -> String {
        format!(
            "Cayley table for D = {} with h(D) = {}",
            self.discriminant().value(),
            self.class_number()
        )
    }

    fn describe(&self) -> String {
        describe_cayley_table(self)
    }
}

impl Visualizable for QuadraticClassGroupGeneratedSubgroup {
    fn format_compact(&self) -> String {
        format!(
            "⟨{}⟩ in Cl(D = {}): order {}/{}",
            self.generator().format_compact(),
            self.discriminant().value(),
            self.order(),
            self.class_number()
        )
    }

    fn describe(&self) -> String {
        describe_generated_subgroup(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::endomorphisms::{
        binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
        candidate_sets::EndomorphismRingCandidateSet,
        quadratic_ideals::QuadraticPrimeBehavior,
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

    #[test]
    fn binary_quadratic_form_description_reports_coefficients_and_invariants() {
        let form = BinaryQuadraticForm::new(1.into(), 1.into(), 6.into());
        let negative_middle = BinaryQuadraticForm::new(2.into(), (-1).into(), 3.into());

        let description = describe_binary_quadratic_form(&form);

        assert!(description.contains("form: (1,1,6)"));
        assert!(description.contains("polynomial: 1x² + 1xy + 6y²"));
        assert!(description.contains("discriminant Δ = b² − 4ac: -23"));
        assert!(description.contains("primitive: yes"));
        assert!(description.contains("reduced positive definite: yes"));
        assert_eq!(form.format_compact(), "(1,1,6)");
        assert!(
            negative_middle
                .describe()
                .contains("polynomial: 2x² - 1xy + 3y²")
        );
    }

    #[test]
    fn cayley_table_description_formats_a_non_cyclic_class_group() {
        let table = QuadraticClassGroup::new(QuadraticDiscriminant::new(-84))
            .expect("D = -84 should be supported")
            .cayley_table()
            .expect("D = -84 should have an enumerated Cayley table");

        let text = describe_cayley_table(&table);

        assert!(text.contains("Quadratic class-group Cayley table"));
        assert!(text.contains("discriminant D: -84"));
        assert!(text.contains("class number h(D): 4"));
        assert!(text.contains("e = (1,0,21)"));
        assert!(text.contains("ab = (5,4,5)"));
        assert!(text.contains("  * |  e  a  b ab"));
        assert!(text.contains("  a |  a  e ab  b"));
        assert!(text.contains("construction cost: Θ(h(D)²) class compositions"));
        assert_eq!(
            table.format_compact(),
            "Cayley table for D = -84 with h(D) = 4"
        );
    }

    #[test]
    fn cayley_table_description_uses_cyclic_power_labels_when_possible() {
        let table = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
            .expect("D = -23 should be supported")
            .cayley_table()
            .expect("D = -23 should have an enumerated Cayley table");

        let text = table.describe();

        assert!(text.contains("g = (2,-1,3)"));
        assert!(text.contains("g² = (2,1,3)"));
        assert!(text.contains("  g |  g g²  e"));
    }

    #[test]
    fn generated_subgroup_description_marks_whether_it_is_the_whole_group() {
        let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
            .expect("D = -23 should be supported");
        let subgroup = class_group
            .generated_subgroup(&BinaryQuadraticForm::new(2.into(), (-1).into(), 3.into()))
            .expect("D = -23 generator should produce a subgroup");

        let text = subgroup.describe();

        assert!(text.contains("Quadratic class-group generated subgroup"));
        assert!(text.contains("generator g: (2,-1,3)"));
        assert!(text.contains("subgroup order: 3"));
        assert!(text.contains("ambient class number h(D): 3"));
        assert!(text.contains("generates whole class group: yes"));
        assert!(text.contains("g⁰ = (1,1,6)"));
        assert!(text.contains("g² = (2,1,3)"));
        assert!(text.contains("not a certified crater action"));
        assert_eq!(
            subgroup.format_compact(),
            "⟨(2,-1,3)⟩ in Cl(D = -23): order 3/3"
        );
    }

    #[test]
    fn generated_subgroup_description_handles_proper_subgroups() {
        let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-84))
            .expect("D = -84 should be supported");
        let subgroup = class_group
            .generated_subgroup(&BinaryQuadraticForm::new(2.into(), 2.into(), 11.into()))
            .expect("D = -84 order-two element should produce a subgroup");

        let text = subgroup.describe();

        assert!(text.contains("subgroup order: 2"));
        assert!(text.contains("ambient class number h(D): 4"));
        assert!(text.contains("generates whole class group: no"));
    }
}
