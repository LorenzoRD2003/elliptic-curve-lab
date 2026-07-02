use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::*;

type F43 = crate::fields::Fp43;

crate::fields::extension_field::define_fp_quadratic_extension!(
    spec: F43Sqrt2Spec,
    field: F43Sqrt2,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2))",
);

#[test]
fn torsion_reports_capture_relative_and_absolute_behavior() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid base curve");
    let relative = base_curve
        .relative_frobenius_on_exact_torsion(2)
        .expect("relative Frobenius on exact two-torsion should evaluate");
    assert!(relative.all_fixed());

    let lifted_curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::zero()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let absolute = lifted_curve
        .absolute_frobenius_on_exact_torsion(4, 1)
        .expect("absolute Frobenius on exact four-torsion should evaluate");
    assert!(
        absolute
            .points()
            .iter()
            .any(|point| !point.fixed_by_frobenius())
    );
}
