use crate::elliptic_curves::{ShortWeierstrassCurve, traits::HasJInvariant};
use crate::fields::{Q, traits::Field};

use super::shared::{F7, f7_curve, q};

#[test]
fn accessors_discriminant_and_rhs_match_the_model() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(F7::eq(curve.a(), &F7::from_i64(2)));
    assert!(F7::eq(curve.b(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c4(), &F7::from_i64(2)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.rhs_value(&F7::from_i64(2)), &F7::from_i64(1)));
}

#[test]
fn weierstrass_invariants_satisfy_the_classical_relation() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");
    let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
    let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

    assert!(F7::eq(&left, &right));
}

#[test]
fn j_invariant_matches_a_classical_exact_example_over_q() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

    assert!(Q::eq(&curve.discriminant(), &q(64, 1)));
    assert!(Q::eq(&curve.c4(), &q(48, 1)));
    assert!(Q::eq(&curve.c6(), &q(0, 1)));
    assert!(Q::eq(&curve.j_invariant(), &q(1728, 1)));
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_method() {
    let curve = f7_curve();

    assert!(F7::eq(
        &HasJInvariant::j_invariant(&curve),
        &ShortWeierstrassCurve::j_invariant(&curve),
    ));
}
