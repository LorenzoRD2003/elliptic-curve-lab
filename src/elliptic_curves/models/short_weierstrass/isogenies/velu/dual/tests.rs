use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isogenies::{DualVeluIsogeny, VeluIsogeny},
    traits::{AffineCurveModel, EnumerableCurveModel},
};
use crate::isogenies::{scalar_multiplication::ScalarMultiplicationIsogeny, traits::Isogeny};

type F41 = crate::fields::Fp41;
type F29 = crate::fields::Fp29;
type Curve = ShortWeierstrassCurve<F41>;
type CurveF29 = ShortWeierstrassCurve<F29>;
type Dual = DualVeluIsogeny<F41>;

fn curve() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn sample_phi() -> VeluIsogeny<Curve> {
    let curve = curve();
    let generator = curve
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample generator should lie on the curve");
    VeluIsogeny::from_generator(curve, generator).expect("sample Vélu isogeny should build")
}

fn curve_f29() -> CurveF29 {
    CurveF29::new(F29::from_i64(2), F29::from_i64(2)).expect("valid F29 curve")
}

fn sample_degree_three_phi() -> VeluIsogeny<CurveF29> {
    let curve = curve_f29();
    let generator = curve
        .point(F29::from_i64(10), F29::from_i64(23))
        .expect("sample degree-three generator should lie on the curve");
    VeluIsogeny::from_generator(curve, generator)
        .expect("sample degree-three Vélu isogeny should build")
}

fn left_relation_holds(phi: &VeluIsogeny<Curve>, dual: &Dual) -> bool {
    let scalar =
        ScalarMultiplicationIsogeny::new(phi.domain().clone(), phi.degree() as u64).unwrap();
    phi.domain().points().into_iter().all(|point| {
        dual.evaluate(&phi.evaluate(&point).unwrap()).unwrap() == scalar.evaluate(&point).unwrap()
    })
}

fn right_relation_holds(phi: &VeluIsogeny<Curve>, dual: &Dual) -> bool {
    let scalar =
        ScalarMultiplicationIsogeny::new(phi.codomain().clone(), phi.degree() as u64).unwrap();
    phi.codomain().points().into_iter().all(|point| {
        phi.evaluate(&dual.evaluate(&point).unwrap()).unwrap() == scalar.evaluate(&point).unwrap()
    })
}

#[test]
fn dual_search_finds_a_degree_matching_candidate_on_the_f41_example() {
    let phi = sample_phi();
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found on the sample curve");
    assert_eq!(dual.degree(), phi.degree());
}

#[test]
fn dual_search_candidate_satisfies_both_duality_relations_on_rational_points() {
    let phi = sample_phi();
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found on the sample curve");
    assert!(left_relation_holds(&phi, &dual));
    assert!(right_relation_holds(&phi, &dual));
}

#[test]
fn public_left_dual_relation_verifier_accepts_the_enumerated_dual() {
    let phi = sample_phi();
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found on the sample curve");
    dual.verify_left_dual_relation(&phi)
        .expect("the exhaustively found dual should satisfy the left relation");
}

#[test]
fn public_right_dual_relation_verifier_accepts_the_enumerated_dual() {
    let phi = sample_phi();
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found on the sample curve");
    dual.verify_right_dual_relation(&phi)
        .expect("the exhaustively found dual should satisfy the right relation");
}

#[test]
fn degree_three_dual_search_reports_an_honest_outcome() {
    let phi = sample_degree_three_phi();
    if let Ok(dual) = phi.find_dual_exhaustively() {
        assert_eq!(dual.degree(), phi.degree());
        dual.verify_left_dual_relation(&phi)
            .expect("a reported dual should satisfy the left relation");
        dual.verify_right_dual_relation(&phi)
            .expect("a reported dual should satisfy the right relation");
    }
}
