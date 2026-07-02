use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::torsion::matrix::{FrobeniusTorsionMatrixError, NTorsionBasis},
    traits::{AffineCurveModel, FiniteGroupCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::traits::*;

type F17 = crate::fields::Fp17;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F17Sqrt3Spec,
    field: F17Sqrt3,
    base: F17,
    non_residue: 3,
    name: "F17(√3)",
);

fn f5_noncyclic_curve() -> ShortWeierstrassCurve<crate::fields::Fp5> {
    ShortWeierstrassCurve::<crate::fields::Fp5>::new(
        crate::fields::Fp5::from_i64(-1),
        crate::fields::Fp5::zero(),
    )
    .expect("valid F5 curve")
}

fn alpha() -> <F17Sqrt3 as Field>::Elem {
    F17Sqrt3::element(vec![F17::zero(), F17::one()])
}

fn lift_f17_curve_to_f17_sqrt3(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<F17Sqrt3> {
    ShortWeierstrassCurve::<F17Sqrt3>::new(
        F17Sqrt3::from_base(*curve.a()),
        F17Sqrt3::from_base(*curve.b()),
    )
    .expect("lifting an F17 curve to F17² should preserve smoothness")
}

fn find_f17_curve_with_nontrivial_two_torsion_frobenius_basis() -> (
    ShortWeierstrassCurve<F17>,
    ShortWeierstrassCurve<F17Sqrt3>,
    NTorsionBasis<crate::elliptic_curves::AffinePoint<F17Sqrt3>>,
) {
    let base_curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(-3), F17::zero())
        .expect("y² = x³ - 3x should be smooth over F17");
    let lifted_curve = lift_f17_curve_to_f17_sqrt3(&base_curve);
    let zero_point = lifted_curve
        .point(F17Sqrt3::zero(), F17Sqrt3::zero())
        .expect("(0,0) should be 2-torsion");
    let alpha_point = lifted_curve
        .point(alpha(), F17Sqrt3::zero())
        .expect("(√3,0) should be 2-torsion");
    let basis = NTorsionBasis::new(&lifted_curve, 2, zero_point, alpha_point)
        .expect("two distinct nonzero 2-torsion points should form a basis");
    (base_curve, lifted_curve, basis)
}

#[test]
fn rational_two_torsion_basis_over_the_base_field_gives_the_identity_matrix() {
    let curve = f5_noncyclic_curve();
    let two_torsion = curve
        .points_of_exact_order(2)
        .expect("F5 curve should have full rational 2-torsion");
    let basis = NTorsionBasis::new(&curve, 2, two_torsion[0].clone(), two_torsion[1].clone())
        .expect("two independent rational 2-torsion points should form a basis");
    let trace = curve
        .frobenius_trace()
        .expect("small enumerable F5 curve should supply a Frobenius trace");
    let report = curve
        .frobenius_matrix_on_n_torsion_basis(trace, basis)
        .expect("matrix report should build on rational 2-torsion");

    assert_eq!(report.matrix().entries(), [[1, 0], [0, 1]]);
}

#[test]
fn nontrivial_extension_two_torsion_basis_still_matches_trace_and_degree_mod_n() {
    let (base_curve, lifted_curve, basis) =
        find_f17_curve_with_nontrivial_two_torsion_frobenius_basis();
    let trace = base_curve
        .frobenius_trace()
        .expect("base F17 curve should supply a Frobenius trace");
    let report = lifted_curve
        .frobenius_matrix_on_n_torsion_basis(trace.clone(), basis)
        .expect("matrix report should build over the lifted curve");

    assert_eq!(report.matrix().modulus(), 2);
    assert!(report.trace_matches_mod_n());
    assert!(report.determinant_matches_mod_n());
}

#[test]
fn torsion_basis_rejects_dependent_points() {
    let curve = f5_noncyclic_curve();
    let two_torsion = curve
        .points_of_exact_order(2)
        .expect("F5 curve should have full rational 2-torsion");

    assert_eq!(
        NTorsionBasis::new(&curve, 2, two_torsion[0].clone(), two_torsion[0].clone()),
        Err(FrobeniusTorsionMatrixError::DependentTorsionBasis)
    );
}
