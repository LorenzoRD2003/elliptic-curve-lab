use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::*;

type F43 = crate::fields::Fp43;
type F17 = crate::fields::Fp17;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F17Sqrt3Spec,
    field: F17Sqrt3,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3))",
);

#[test]
fn torsion_reports_capture_relative_and_absolute_behavior() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid base curve");
    let relative = base_curve
        .relative_frobenius_on_exact_torsion(2)
        .expect("relative Frobenius on exact two-torsion should evaluate");
    assert!(relative.all_fixed());

    let lifted_curve = ShortWeierstrassCurve::<F17Sqrt3>::new(
        F17Sqrt3::from_base(F17::zero()),
        F17Sqrt3::from_base(F17::one()),
    )
    .expect("base-defined curve should stay smooth over F17^2");
    let absolute = lifted_curve
        .absolute_frobenius_on_exact_torsion(2, 1)
        .expect("absolute Frobenius on exact two-torsion should evaluate");
    assert!(
        absolute
            .points()
            .iter()
            .any(|point| !point.fixed_by_frobenius())
    );
}
