use crate::isogenies::class_group_action::{
    ClassGroupIsogenyActionReport, ClassGroupIsogenyActionSegment, CraterDirectionCertification,
};
use crate::visualization::Visualizable;

impl Visualizable for ClassGroupIsogenyActionReport {
    fn format_compact(&self) -> String {
        format!(
            "class-group isogeny action: v{} -> v{} through {} local segment(s)",
            self.start().0,
            self.target().0,
            self.segments().len()
        )
    }

    fn describe(&self) -> String {
        describe_class_group_isogeny_action(self)
    }
}

fn describe_class_group_isogeny_action(report: &ClassGroupIsogenyActionReport) -> String {
    let mut lines = vec![
        "Class-group isogeny action execution".to_string(),
        "-------------------------------------".to_string(),
        format!("start: v{}", report.start().0),
        format!("target: v{}", report.target().0),
        format!("local segments: {}", report.segments().len()),
        "segments:".to_string(),
    ];

    if report.segments().is_empty() {
        lines.push("  none (principal class)".to_string());
    } else {
        for segment in report.segments() {
            lines.extend(format_segment(segment));
        }
    }

    lines.push(
        "interpretation: this concatenates supplied local oriented crater witnesses; it does not infer arithmetic orientation."
            .to_string(),
    );
    lines.join("\n")
}

fn format_segment(segment: &ClassGroupIsogenyActionSegment) -> Vec<String> {
    vec![
        format!("  factor {}:", segment.factor_index() + 1),
        format!("    ideal norm: {}", segment.ideal().norm()),
        format!(
            "    form class: {}",
            segment.generator_form().format_compact()
        ),
        format!("    exponent: {}", segment.exponent()),
        format!("    path: {}", format_path(segment)),
        format!(
            "    direction: {}",
            format_direction(segment.direction_certification())
        ),
    ]
}

fn format_path(segment: &ClassGroupIsogenyActionSegment) -> String {
    segment
        .path()
        .iter()
        .map(|node| format!("v{}", node.0))
        .collect::<Vec<_>>()
        .join(" -> ")
}

fn format_direction(certification: CraterDirectionCertification) -> &'static str {
    match certification {
        CraterDirectionCertification::GraphDeterministic => "graph-deterministic",
        CraterDirectionCertification::UserSuppliedArithmeticOrientation => {
            "user-supplied crater orientation"
        }
        CraterDirectionCertification::CertifiedArithmeticOrientation => "certified arithmetic",
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use num_bigint::BigInt;

    use super::*;
    use crate::elliptic_curves::endomorphisms::{
        binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    };
    use crate::fields::Fp101;
    use crate::isogenies::{
        class_group_action::{ClassGroupActionPlan, CraterOrientationWitness},
        graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
    };
    use crate::{elliptic_curves::ShortWeierstrassCurve, visualization::Visualizable};

    #[test]
    fn executed_action_description_stays_witness_based() {
        let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
            .expect("D = -23 should define a class group");
        let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), 1u8.into())
            .expect("D = -23 should define an imaginary quadratic order");
        let ideal = PrimeNormIdeal::split(order, 3u8.into(), 1u8.into()).expect("3 should split");
        let target = BinaryQuadraticForm::new(BigInt::from(2), BigInt::from(-1), BigInt::from(3));
        let plan = ClassGroupActionPlan::from_local_ideals(
            &class_group,
            std::slice::from_ref(&ideal),
            &target,
        )
        .expect("local generator should produce a plan");
        let curve = ShortWeierstrassCurve::<Fp101>::new(Fp101::from_i64(1), Fp101::from_i64(12))
            .expect("valid F_101 curve");
        let graph = IsogenyGraphBuilder::new(curve, 3)
            .max_depth(3)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small graph should build");
        let crater = graph
            .volcano_crater_report(ideal.norm())
            .expect("crater report should build");
        let labeled = graph
            .labeled_crater_walk_report(&class_group, ideal, IsogenyGraphNodeId(0))
            .expect("labeled walk should build");
        let orientation = CraterOrientationWitness::new(
            &crater,
            BTreeMap::from([
                (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
                (IsogenyGraphNodeId(1), IsogenyGraphNodeId(2)),
                (IsogenyGraphNodeId(2), IsogenyGraphNodeId(0)),
            ]),
        )
        .expect("crater orientation should validate");
        let oriented = labeled
            .with_user_orientation(orientation)
            .expect("orientation should attach");

        let report = plan
            .execute_from(IsogenyGraphNodeId(0), &[oriented])
            .expect("action should execute");
        let text = report.describe();

        assert!(text.contains("Class-group isogeny action execution"));
        assert!(text.contains("path: v0 -> v1"));
        assert!(text.contains("user-supplied crater orientation"));
        assert!(text.contains("does not infer arithmetic orientation"));
    }
}
