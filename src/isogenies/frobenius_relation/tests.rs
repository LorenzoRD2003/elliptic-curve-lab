use proptest::prelude::*;

use crate::elliptic_curves::short_weierstrass::isogenies::VeluIsogeny;
use crate::elliptic_curves::{ShortWeierstrassCurve, traits::AffineCurveModel};
use crate::fields::{Fp, traits::Field};
use crate::isogenies::{
    frobenius_relation::{FrobeniusComparableIsogeny, FrobeniusComparableIsogenyGraph},
    graphs::IsogenyGraphBuilder,
    scalar_multiplication::ScalarMultiplicationIsogeny,
};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};

type F41 = Fp<41>;

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

#[test]
fn scalar_and_graph_isogeny_frobenius_checks_hold_on_small_examples() {
    let curve = f41_curve();
    let scalar_isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 2)
        .expect("small scalar-multiplication isogeny should build");
    let scalar_relation = scalar_isogeny
        .frobenius_relation_report()
        .expect("scalar Frobenius relation should compute");
    assert!(scalar_relation.holds());

    let generator = curve
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample generator should lie on the curve");
    let velu =
        VeluIsogeny::from_generator(curve.clone(), generator).expect("Velu isogeny should build");
    let velu_relation = velu
        .frobenius_relation_report()
        .expect("Velu Frobenius relation should compute");
    assert!(velu_relation.holds());

    let graph = IsogenyGraphBuilder::new(curve, 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");
    let graph_report = graph
        .frobenius_relation_report()
        .expect("graph Frobenius relation should compute");
    assert!(graph_report.holds());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[test]
    fn property_scalar_isogenies_preserve_curve_order_and_trace(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
        scalar in 1u64..=4,
    ) {
        let isogeny = ScalarMultiplicationIsogeny::new(curve, scalar)
            .expect("small scalar-multiplication isogenies should build");
        let relation = isogeny
            .frobenius_relation_report()
            .expect("Frobenius relation should compute for scalar isogenies");

        prop_assert!(relation.same_curve_order());
        prop_assert!(relation.same_trace());
        prop_assert!(relation.holds());
    }
}
