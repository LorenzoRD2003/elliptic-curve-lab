use crate::elliptic_curves::endomorphisms::quadratic_ideals::IdealFormConvention;
use crate::isogenies::class_group_action::{
    CraterDirectionCertification, CraterIdealPrimeBehavior, CraterWalkReport,
    CraterWalkTermination, LabeledCraterWalkReport, OrientedCraterClassOrderComparison,
    OrientedCraterClassOrderStatus, OrientedCraterPowerActionReport,
    OrientedLabeledCraterWalkReport,
};
use crate::visualization::{Visualizable, shared::yes_no};

impl Visualizable for LabeledCraterWalkReport {
    fn format_compact(&self) -> String {
        format!(
            "labeled crater walk: {} via {}",
            self.walk().format_compact(),
            self.form_label().reduced_form().format_compact()
        )
    }

    fn describe(&self) -> String {
        describe_labeled_crater_walk(self)
    }
}

impl Visualizable for OrientedLabeledCraterWalkReport {
    fn format_compact(&self) -> String {
        format!("user-oriented {}", self.labeled_walk().format_compact())
    }

    fn describe(&self) -> String {
        describe_oriented_labeled_crater_walk(self)
    }
}

impl Visualizable for OrientedCraterPowerActionReport {
    fn format_compact(&self) -> String {
        format!(
            "oriented crater power n = {}: {}",
            self.exponent(),
            format_nodes(self.path())
        )
    }

    fn describe(&self) -> String {
        describe_oriented_crater_power(self)
    }
}

impl Visualizable for OrientedCraterClassOrderComparison {
    fn format_compact(&self) -> String {
        format!(
            "class order {} vs oriented orbit {} ({})",
            self.class_order(),
            self.oriented_orbit_length()
                .map(|length| length.to_string())
                .unwrap_or_else(|| "not closed".to_string()),
            format_class_order_status(self.status())
        )
    }

    fn describe(&self) -> String {
        describe_class_order_comparison(self)
    }
}

fn describe_labeled_crater_walk(report: &LabeledCraterWalkReport) -> String {
    let local_label = report.local_label();
    let form_label = report.form_label();
    let walk = report.walk();

    let mut lines = vec![
        "Labeled crater walk".to_string(),
        "-------------------".to_string(),
        format!("local ideal label: {}", format_prime_ideal_label(report)),
        format!("local prime ℓ: {}", local_label.crater_prime()),
        format!(
            "prime behavior: {}",
            format_prime_behavior(local_label.prime_behavior())
        ),
        format!(
            "associated raw form: {}",
            form_label.raw_form().format_compact()
        ),
        format!(
            "associated form class: {}",
            form_label.reduced_form().format_compact()
        ),
        format!(
            "form convention: {}",
            format_form_convention(form_label.convention())
        ),
        format!(
            "direction: {}",
            format_direction_certification(report.direction_certification())
        ),
        "arithmetic orientation: not certified as arithmetic orientation".to_string(),
        format!(
            "conjugate direction: {}",
            format_conjugate_direction(local_label.prime_behavior())
        ),
        format!("visited: {}", format_node_path(walk)),
        format!(
            "cycle length: {}",
            walk.cycle_length()
                .map(|length| length.to_string())
                .unwrap_or_else(|| "not closed".to_string())
        ),
        format!(
            "start in certified crater: {}",
            yes_no(walk.start_in_crater())
        ),
        format!(
            "graph termination: {}",
            format_walk_termination(walk.termination())
        ),
    ];

    lines.push(
        "Interpretation: the ideal and form labels are compatible with the crater; this is not an arithmetic action computation."
            .to_string(),
    );

    lines.join("\n")
}

fn describe_oriented_labeled_crater_walk(report: &OrientedLabeledCraterWalkReport) -> String {
    let labeled = report.labeled_walk();
    let start = labeled.walk().start();
    let oriented_path = report
        .orientation()
        .oriented_cycle_from(start)
        .map(|nodes| format_nodes(&nodes))
        .unwrap_or_else(|| "not available".to_string());

    [
        "User-supplied crater orientation".to_string(),
        "--------------------------------".to_string(),
        format!(
            "base labeled walk: {}",
            labeled.format_compact()
        ),
        format!(
            "direction: {}",
            format_direction_certification(report.direction_certification())
        ),
        format!("oriented cycle from start: {oriented_path}"),
        "orientation source: user-supplied witness checked against certified crater edges"
            .to_string(),
        "Interpretation: this records a declared crater direction; it is not an automatically certified arithmetic orientation and still does not compute an ideal action on the curve."
            .to_string(),
    ]
    .join("\n")
}

fn describe_oriented_crater_power(report: &OrientedCraterPowerActionReport) -> String {
    [
        "Oriented crater local power".to_string(),
        "---------------------------".to_string(),
        format!("exponent n: {}", report.exponent()),
        format!("generator ideal norm: {}", report.generator_ideal().norm()),
        format!(
            "generator form class: {}",
            report.generator_form().format_compact()
        ),
        format!("oriented path: {}", format_nodes(report.path())),
        format!("target: v{}", report.target().0),
        "modeled local power under user-supplied orientation".to_string(),
        "Interpretation: this walks inside the oriented crater; it does not compare with algebraic class orders."
            .to_string(),
    ]
    .join("\n")
}

fn describe_class_order_comparison(report: &OrientedCraterClassOrderComparison) -> String {
    [
        "Crater class-order comparison".to_string(),
        "-----------------------------".to_string(),
        format!("start: v{}", report.start().0),
        format!("generator ideal norm: {}", report.generator_ideal().norm()),
        format!(
            "generator form class: {}",
            report.generator_form().format_compact()
        ),
        format!("class-group order: {}", report.class_order()),
        format!(
            "oriented orbit length: {}",
            report
                .oriented_orbit_length()
                .map(|length| length.to_string())
                .unwrap_or_else(|| "not closed".to_string())
        ),
        format!("status: {}", format_class_order_status(report.status())),
        "orientation source: user-supplied witness checked against certified crater edges"
            .to_string(),
        "Interpretation: this is a diagnostic comparison; equality would not by itself certify an arithmetic CM action, and a difference is expected when the toy crater does not model the full class orbit."
            .to_string(),
    ]
    .join("\n")
}

fn format_prime_ideal_label(report: &LabeledCraterWalkReport) -> String {
    let ideal = report.local_label().ideal();
    format!("𝔭 = ({}, ω - {})", ideal.norm(), ideal.root_mod_ell())
}

fn format_prime_behavior(behavior: CraterIdealPrimeBehavior) -> &'static str {
    match behavior {
        CraterIdealPrimeBehavior::Split => "split",
        CraterIdealPrimeBehavior::Ramified => "ramified",
    }
}

fn format_conjugate_direction(behavior: CraterIdealPrimeBehavior) -> &'static str {
    match behavior {
        CraterIdealPrimeBehavior::Split => "not distinguished",
        CraterIdealPrimeBehavior::Ramified => "same ramified ideal",
    }
}

fn format_direction_certification(certification: CraterDirectionCertification) -> &'static str {
    match certification {
        CraterDirectionCertification::GraphDeterministic => "graph-deterministic",
        CraterDirectionCertification::UserSuppliedArithmeticOrientation => {
            "user-supplied crater orientation"
        }
        CraterDirectionCertification::CertifiedArithmeticOrientation => {
            "certified arithmetic orientation"
        }
    }
}

fn format_form_convention(convention: IdealFormConvention) -> &'static str {
    match convention {
        IdealFormConvention::RepresentsIdeal => "represents the supplied ideal class",
        IdealFormConvention::RepresentsInverseIdeal => "represents the inverse ideal class",
    }
}

fn format_walk_termination(termination: CraterWalkTermination) -> &'static str {
    match termination {
        CraterWalkTermination::ClosedCycle => "closed cycle",
        CraterWalkTermination::StartOutsideCrater => "start outside certified crater",
        CraterWalkTermination::NoCertifiedOutgoingEdge => "no certified outgoing crater edge",
        CraterWalkTermination::RepeatedNonStartNode => "repeated non-start crater node",
    }
}

fn format_class_order_status(status: OrientedCraterClassOrderStatus) -> &'static str {
    match status {
        OrientedCraterClassOrderStatus::MatchesOrientedOrbit => "matches oriented orbit",
        OrientedCraterClassOrderStatus::OrientedOrbitLengthDiffers => {
            "oriented orbit length differs"
        }
        OrientedCraterClassOrderStatus::OrbitDidNotClose => "oriented orbit did not close",
    }
}

fn format_node_path(report: &CraterWalkReport) -> String {
    format_nodes(report.visited())
}

fn format_nodes(nodes: &[crate::isogenies::graphs::IsogenyGraphNodeId]) -> String {
    nodes
        .iter()
        .map(|node| format!("v{}", node.0))
        .collect::<Vec<_>>()
        .join(" -> ")
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        endomorphisms::{
            binary_quadratic_forms::QuadraticClassGroup,
            quadratic_ideals::PrimeNormIdeal,
            quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
        },
    };
    use crate::fields::Fp7;
    use crate::isogenies::{
        class_group_action::{
            CraterOrientationWitness, LabeledCraterWalkReport, OrientedCraterClassOrderStatus,
        },
        graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
    };
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

    fn class_group_minus_23() -> QuadraticClassGroup {
        QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
            .expect("D = -23 should define an imaginary quadratic class group")
    }

    #[test]
    fn labeled_crater_walk_explanation_keeps_graph_and_arithmetic_claims_separate() {
        let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let ideal = split_three_ideal();
        let crater = graph
            .volcano_crater_report(ideal.norm())
            .expect("crater report should build for the ideal norm");
        let report = LabeledCraterWalkReport::from_crater_report(
            &crater,
            &class_group_minus_23(),
            ideal,
            IsogenyGraphNodeId(0),
        )
        .expect("matching crater, ideal, and class group should label the walk");

        let explanation = report.describe();

        assert!(explanation.contains("Labeled crater walk"));
        assert!(explanation.contains("local ideal label"));
        assert!(explanation.contains("associated form class"));
        assert!(explanation.contains("(2,-1,3)"));
        assert!(explanation.contains("prime behavior: split"));
        assert!(explanation.contains("graph-deterministic"));
        assert!(explanation.contains("not certified as arithmetic orientation"));
        assert!(explanation.contains("visited: v0 -> v1 -> v0"));
        assert!(explanation.contains("cycle length: 2"));
        assert!(explanation.contains("graph termination"));
        assert!(!explanation.contains("class order"));
        assert!(!explanation.contains("computed action"));
        assert!(!explanation.contains("computed class-group action"));
        assert!(!explanation.contains("[𝔭] * E"));
        assert!(report.format_compact().contains("labeled crater walk"));
    }

    #[test]
    fn oriented_labeled_crater_walk_explanation_mentions_user_supplied_orientation() {
        let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let ideal = split_three_ideal();
        let crater = graph
            .volcano_crater_report(ideal.norm())
            .expect("crater report should build for the ideal norm");
        let report = LabeledCraterWalkReport::from_crater_report(
            &crater,
            &class_group_minus_23(),
            ideal,
            IsogenyGraphNodeId(0),
        )
        .expect("matching crater, ideal, and class group should label the walk");
        let witness = CraterOrientationWitness::new(
            &crater,
            [
                (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
                (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
            ]
            .into(),
        )
        .expect("two-node crater has a certified orientation witness");

        let oriented = report
            .with_user_orientation(witness)
            .expect("witness should attach to labeled walk");
        let explanation = oriented.describe();

        assert!(explanation.contains("User-supplied crater orientation"));
        assert!(explanation.contains("user-supplied crater orientation"));
        assert!(explanation.contains("oriented cycle from start: v0 -> v1 -> v0"));
        assert!(explanation.contains("user-supplied witness"));
        assert!(!explanation.contains("computed class-group action"));
        assert!(!explanation.contains("[𝔭] * E"));
    }

    #[test]
    fn oriented_crater_power_explanation_shows_exponent_path_and_target() {
        let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let ideal = split_three_ideal();
        let crater = graph
            .volcano_crater_report(ideal.norm())
            .expect("crater report should build for the ideal norm");
        let report = LabeledCraterWalkReport::from_crater_report(
            &crater,
            &class_group_minus_23(),
            ideal,
            IsogenyGraphNodeId(0),
        )
        .expect("matching crater, ideal, and class group should label the walk");
        let witness = CraterOrientationWitness::new(
            &crater,
            [
                (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
                (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
            ]
            .into(),
        )
        .expect("two-node crater has a certified orientation witness");
        let oriented = report
            .with_user_orientation(witness)
            .expect("witness should attach to labeled walk");

        let power = oriented
            .apply_power_from(IsogenyGraphNodeId(0), 1.into())
            .expect("local power should apply");
        let explanation = power.describe();

        assert!(explanation.contains("Oriented crater local power"));
        assert!(explanation.contains("exponent n: 1"));
        assert!(explanation.contains("generator form class: (2,-1,3)"));
        assert!(explanation.contains("oriented path: v0 -> v1"));
        assert!(explanation.contains("target: v1"));
        assert!(explanation.contains("modeled local power under user-supplied orientation"));
        assert!(!explanation.contains("class-group order"));
    }

    #[test]
    fn class_order_comparison_explanation_stays_diagnostic() {
        let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let ideal = split_three_ideal();
        let crater = graph
            .volcano_crater_report(ideal.norm())
            .expect("crater report should build for the ideal norm");
        let report = LabeledCraterWalkReport::from_crater_report(
            &crater,
            &class_group_minus_23(),
            ideal,
            IsogenyGraphNodeId(0),
        )
        .expect("matching crater, ideal, and class group should label the walk");
        let witness = CraterOrientationWitness::new(
            &crater,
            [
                (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
                (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
            ]
            .into(),
        )
        .expect("two-node crater has a certified orientation witness");
        let oriented = report
            .with_user_orientation(witness)
            .expect("witness should attach to labeled walk");

        let comparison = oriented
            .compare_generator_order(&class_group_minus_23(), IsogenyGraphNodeId(0))
            .expect("comparison should compute");
        let explanation = comparison.describe();

        assert_eq!(
            comparison.status(),
            OrientedCraterClassOrderStatus::OrientedOrbitLengthDiffers
        );
        assert!(explanation.contains("Crater class-order comparison"));
        assert!(explanation.contains("generator form class: (2,-1,3)"));
        assert!(explanation.contains("class-group order: 3"));
        assert!(explanation.contains("oriented orbit length: 2"));
        assert!(explanation.contains("oriented orbit length differs"));
        assert!(explanation.contains("diagnostic comparison"));
        assert!(!explanation.contains("certifies an arithmetic CM action"));
        assert!(comparison.format_compact().contains("class order 3"));
    }
}
