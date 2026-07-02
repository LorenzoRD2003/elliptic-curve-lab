use crate::elliptic_curves::short_weierstrass::isogenies::{
    DualVeluIsogeny, VeluIsogeny,
    frobenius::{AbsoluteFrobeniusIsogeny, VerschiebungCertificate, VerschiebungIsogeny},
    function_field_maps::ShortWeierstrassFunctionFieldMap,
};
use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::isogenies::frobenius::FrobeniusLikeIsogeny,
    traits::AffineCurveModel,
};
use crate::isogenies::{
    dual_report::{DualIsogenyReport, DualityKind},
    traits::Isogeny,
};
use num_bigint::BigUint;

type F29 = crate::fields::Fp29;
type F41 = crate::fields::Fp41;
type Curve29 = ShortWeierstrassCurve<F29>;
type Curve41 = ShortWeierstrassCurve<F41>;

fn bu(value: usize) -> BigUint {
    BigUint::from(value)
}

fn curve_f29() -> Curve29 {
    Curve29::new(F29::from_i64(2), F29::from_i64(2)).expect("valid curve")
}

fn curve_f41() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn degree_three_phi() -> VeluIsogeny<Curve29> {
    let curve = curve_f29();
    let generator = curve
        .point(F29::from_i64(10), F29::from_i64(23))
        .expect("sample generator should lie on the curve");
    VeluIsogeny::from_generator(curve, generator).expect("sample Velu isogeny should build")
}

#[test]
fn classical_dual_report_records_basic_duality_data() {
    let phi = degree_three_phi();
    let dual = phi.find_dual_exhaustively().expect("dual should be found");
    let report: DualIsogenyReport<Curve29, Curve29> =
        DualVeluIsogeny::dual_report(&phi, &dual).expect("report should build");

    assert_eq!(report.duality_kind(), DualityKind::SeparableClassical);
    assert_eq!(report.phi_degree(), &bu(3));
    assert_eq!(report.dual_degree(), &bu(3));
    assert!(report.left_relation_holds());
    assert!(report.right_relation_holds());
    assert_eq!(report.phi_kernel_summary().total_degree(), Some(3));
    assert_eq!(report.dual_kernel_summary().total_degree(), Some(3));
    assert_eq!(
        report.phi_degree_factorization().separable_degree(),
        Some(&bu(3))
    );
    assert_eq!(
        report.phi_degree_factorization().inseparable_degree(),
        Some(&bu(1))
    );
    assert_eq!(
        report.dual_degree_factorization().separable_degree(),
        Some(&bu(3))
    );
    assert_eq!(
        report.dual_degree_factorization().inseparable_degree(),
        Some(&bu(1))
    );
}

#[test]
fn frobenius_verschiebung_report_records_known_frobenius_data_and_unknown_v_data() {
    let curve = curve_f41();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let candidate_v = ShortWeierstrassFunctionFieldMap::new(
        frobenius.codomain().clone(),
        frobenius.domain().clone(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .x(),
        frobenius
            .as_function_field_map()
            .codomain_function_field()
            .y(),
    )
    .expect("identity candidate on the twist should validate");
    let verschiebung =
        VerschiebungIsogeny::new(frobenius.clone(), candidate_v).expect("candidate should build");
    let expected_left = frobenius
        .as_function_field_map()
        .compose(verschiebung.as_function_field_map())
        .expect("left composition should build");
    let expected_right = verschiebung
        .as_function_field_map()
        .compose(&frobenius.as_function_field_map())
        .expect("right composition should build");
    let certificate = VerschiebungCertificate::new(verschiebung, expected_left, expected_right)
        .expect("certificate should build");
    let report = certificate.dual_report().expect("report should build");

    assert_eq!(report.duality_kind(), DualityKind::FrobeniusVerschiebung);
    assert_eq!(report.phi_degree(), &bu(41));
    assert_eq!(report.dual_degree(), &bu(41));
    assert_eq!(
        report.phi_degree_factorization().separable_degree(),
        Some(&bu(1))
    );
    assert_eq!(
        report.phi_degree_factorization().inseparable_degree(),
        Some(&bu(41))
    );
    assert!(
        report
            .dual_degree_factorization()
            .separable_degree()
            .is_none()
    );
    assert!(report.left_relation_holds());
    assert!(report.right_relation_holds());
    assert!(
        report
            .dual_kernel_summary()
            .short_label()
            .contains("not yet modeled")
    );
}
