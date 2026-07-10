use super::{bu, class_group_minus_23, crater_report, f7_curve, form, split_three_ideal};
use crate::{
    elliptic_curves::endomorphisms::{
        binary_quadratic_forms::QuadraticClassGroup, quadratic_ideals::IdealFormConvention,
        quadratic_orders::QuadraticDiscriminant,
    },
    isogenies::{
        class_group_action::{
            CraterDirectionCertification, CraterIdealLabelError, CraterIdealPrimeBehavior,
            LabeledCraterWalkReport,
        },
        graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId, endomorphisms::CraterShape},
    },
};

#[test]
fn labeled_crater_walk_report_combines_walk_local_label_and_form_label() {
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");
    let ideal = split_three_ideal();
    let class_group = class_group_minus_23();
    let crater = graph
        .volcano_crater_report(ideal.norm())
        .expect("crater report should build for the ideal norm");

    let report = LabeledCraterWalkReport::from_crater_report(
        &crater,
        &class_group,
        ideal,
        IsogenyGraphNodeId(0),
    )
    .expect("matching crater, ideal, and class group should produce a labeled walk");

    assert_eq!(report.local_label().crater_prime(), &bu(3));
    assert_eq!(report.local_label().ideal().norm(), &bu(3));
    assert_eq!(
        report.local_label().prime_behavior(),
        CraterIdealPrimeBehavior::Split
    );
    assert_eq!(report.walk().start(), IsogenyGraphNodeId(0));
    assert_eq!(
        report.walk().crater_shape(),
        CraterShape::TwoVertex {
            directed_edge_count: 2
        }
    );
    assert_eq!(
        report.walk().visited(),
        &[
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(0)
        ]
    );
    assert_eq!(report.walk().cycle_length(), Some(2));
    assert_eq!(
        report.form_label().raw_form(),
        &form(3, 1, 2),
        "root 1 over D = -23 gives the GP/PARI raw form (3,1,2)"
    );
    assert_eq!(report.form_label().reduced_form(), &form(2, -1, 3));
    assert!(class_group.contains_reduced_form(report.form_label().reduced_form()));
    assert_eq!(
        report.form_label().convention(),
        IdealFormConvention::RepresentsIdeal
    );
    assert_eq!(
        report.direction_certification(),
        CraterDirectionCertification::GraphDeterministic
    );
}

#[test]
fn labeled_crater_walk_report_rejects_crater_prime_mismatch() {
    let crater = crater_report(bu(5), Vec::new());
    let ideal = split_three_ideal();

    let error = LabeledCraterWalkReport::from_crater_report(
        &crater,
        &class_group_minus_23(),
        ideal,
        IsogenyGraphNodeId(0),
    )
    .expect_err("crater prime should match the ideal norm before labeling the walk");

    assert_eq!(
        error,
        CraterIdealLabelError::PrimeNormMismatch {
            ideal_norm: bu(3),
            crater_prime: bu(5),
        }
    );
}

#[test]
fn labeled_crater_walk_report_rejects_class_group_discriminant_mismatch() {
    let crater = crater_report(bu(3), Vec::new());
    let ideal = split_three_ideal();
    let wrong_class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-20))
        .expect("D = -20 should define an imaginary quadratic class group");

    let error = LabeledCraterWalkReport::from_crater_report(
        &crater,
        &wrong_class_group,
        ideal,
        IsogenyGraphNodeId(0),
    )
    .expect_err("class group discriminant should match before labeling the walk");

    assert_eq!(
        error,
        CraterIdealLabelError::OrderDiscriminantMismatch {
            ideal_discriminant: (-23).into(),
            class_group_discriminant: (-20).into(),
        }
    );
}
