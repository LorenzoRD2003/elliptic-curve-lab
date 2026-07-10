use super::{
    bu, class_group_minus_23, crater_report, ramified_twenty_three_ideal, split_three_ideal,
};
use crate::{
    elliptic_curves::endomorphisms::{
        binary_quadratic_forms::QuadraticClassGroup, quadratic_orders::QuadraticDiscriminant,
    },
    isogenies::class_group_action::{
        CraterIdealLabelError, CraterIdealLabelReport, CraterIdealPrimeBehavior,
    },
};

#[test]
fn crater_ideal_label_accepts_split_ideal_for_matching_crater_and_class_group() {
    let crater = crater_report(bu(3), Vec::new());
    let ideal = split_three_ideal();

    let report = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect("split ideal should label the matching crater");

    assert_eq!(report.crater_prime(), &bu(3));
    assert_eq!(report.ideal().norm(), &bu(3));
    assert_eq!(report.prime_behavior(), CraterIdealPrimeBehavior::Split);
}

#[test]
fn crater_ideal_label_accepts_ramified_ideal_for_matching_crater_and_class_group() {
    let crater = crater_report(bu(23), Vec::new());
    let ideal = ramified_twenty_three_ideal();

    let report = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect("ramified ideal should label the matching crater");

    assert_eq!(report.crater_prime(), &bu(23));
    assert_eq!(report.ideal().norm(), &bu(23));
    assert_eq!(report.prime_behavior(), CraterIdealPrimeBehavior::Ramified);
}

#[test]
fn crater_ideal_label_rejects_crater_prime_mismatch() {
    let crater = crater_report(bu(5), Vec::new());
    let ideal = split_three_ideal();

    let error = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect_err("crater prime should match the ideal norm");

    assert_eq!(
        error,
        CraterIdealLabelError::PrimeNormMismatch {
            ideal_norm: bu(3),
            crater_prime: bu(5),
        }
    );
}

#[test]
fn crater_ideal_label_rejects_class_group_discriminant_mismatch() {
    let crater = crater_report(bu(3), Vec::new());
    let ideal = split_three_ideal();
    let wrong_class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-20))
        .expect("D = -20 should define an imaginary quadratic class group");

    let error = CraterIdealLabelReport::new(&crater, &wrong_class_group, ideal)
        .expect_err("ideal order and class group must have the same discriminant");

    assert_eq!(
        error,
        CraterIdealLabelError::OrderDiscriminantMismatch {
            ideal_discriminant: (-23).into(),
            class_group_discriminant: (-20).into(),
        }
    );
}
