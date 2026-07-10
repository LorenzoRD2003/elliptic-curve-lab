mod class_order_comparison;
mod crater_walk;
mod horizontal_ideal;
mod ideal_label;
mod labeled_crater_walk;
mod orientation;
mod oriented_power;

use num_bigint::BigUint;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::{
        binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    },
};
use crate::isogenies::graphs::{
    IsogenyGraphEdgeId, IsogenyGraphNodeId,
    endomorphisms::{
        CraterReport, CraterShape, HorizontalEdgeReport, HorizontalEdgeStatus,
        VolcanoStructureReport,
    },
};

type F101 = crate::fields::Fp101;

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn form(a: i64, b: i64, c: i64) -> BinaryQuadraticForm {
    BinaryQuadraticForm::new(a.into(), b.into(), c.into())
}

fn order_minus_23() -> ImaginaryQuadraticOrder {
    ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(1))
        .expect("D = -23 should define an imaginary quadratic maximal order")
}

fn cm_field_minus_23_curve() -> ShortWeierstrassCurve<F101> {
    ShortWeierstrassCurve::<F101>::new(F101::from_i64(1), F101::from_i64(12))
        .expect("valid F_101 curve with Frobenius field Q(sqrt(-23))")
}

fn split_three_ideal() -> PrimeNormIdeal {
    PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23")
}

fn ramified_twenty_three_ideal() -> PrimeNormIdeal {
    PrimeNormIdeal::ramified(order_minus_23(), bu(23))
        .expect("23 ramifies in the order of discriminant -23")
}

fn horizontal_edge(status: HorizontalEdgeStatus) -> HorizontalEdgeReport {
    HorizontalEdgeReport::new(
        IsogenyGraphEdgeId(7),
        IsogenyGraphNodeId(1),
        IsogenyGraphNodeId(2),
        status,
    )
}

fn crater_report(prime: BigUint, edges: Vec<HorizontalEdgeReport>) -> CraterReport {
    crater_report_with_nodes(prime, Vec::new(), edges)
}

fn crater_report_with_nodes(
    prime: BigUint,
    crater_nodes: Vec<IsogenyGraphNodeId>,
    edges: Vec<HorizontalEdgeReport>,
) -> CraterReport {
    CraterReport::new(
        prime.clone(),
        VolcanoStructureReport::from_floor_paths(prime, Vec::new(), Vec::new()),
        crater_nodes,
        edges,
        CraterShape::EmptyCertifiedCrater,
    )
}

fn class_group_minus_23() -> QuadraticClassGroup {
    QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should define an imaginary quadratic class group")
}
