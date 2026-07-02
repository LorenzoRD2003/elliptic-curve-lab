use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::torsion::matrix::{FrobeniusTorsionMatrixError, NTorsionBasis},
    traits::{AffineCurveModel, FiniteGroupCurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::traits::*;

type F43 = crate::fields::Fp43;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F43Sqrt2Spec,
    field: F43Sqrt2,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2))",
);

fn f5_noncyclic_curve() -> ShortWeierstrassCurve<crate::fields::Fp5> {
    ShortWeierstrassCurve::<crate::fields::Fp5>::new(
        crate::fields::Fp5::from_i64(-1),
        crate::fields::Fp5::zero(),
    )
    .expect("valid F5 curve")
}

fn alpha() -> <F43Sqrt2 as Field>::Elem {
    F43Sqrt2::element(vec![F43::zero(), F43::one()])
}

fn lift_f43_curve_to_f43_sqrt2(
    curve: &ShortWeierstrassCurve<F43>,
) -> ShortWeierstrassCurve<F43Sqrt2> {
    ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(*curve.a()),
        F43Sqrt2::from_base(*curve.b()),
    )
    .expect("lifting an F43 curve to F43^2 should preserve smoothness")
}

fn find_f43_curve_with_nontrivial_two_torsion_frobenius_basis() -> (
    ShortWeierstrassCurve<F43>,
    ShortWeierstrassCurve<F43Sqrt2>,
    NTorsionBasis<crate::elliptic_curves::AffinePoint<F43Sqrt2>>,
) {
    let base_curve = ShortWeierstrassCurve::<F43>::new(F43::from_i64(-2), F43::zero())
        .expect("y^2 = x^3 - 2x should be smooth over F43");
    let lifted_curve = lift_f43_curve_to_f43_sqrt2(&base_curve);
    let zero_point = lifted_curve
        .point(F43Sqrt2::zero(), F43Sqrt2::zero())
        .expect("(0,0) should be 2-torsion");
    let alpha_point = lifted_curve
        .point(alpha(), F43Sqrt2::zero())
        .expect("(sqrt(2),0) should be 2-torsion");
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
        find_f43_curve_with_nontrivial_two_torsion_frobenius_basis();
    let trace = base_curve
        .frobenius_trace()
        .expect("base F43 curve should supply a Frobenius trace");
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
