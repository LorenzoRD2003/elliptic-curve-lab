use core::fmt;

use crate::elliptic_curves::TwistedEdwardsCurve;
use crate::fields::traits::Field;
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::montgomery::format_montgomery_curve,
    shared::{format_field_elem as format_elem, parenthesize_if_needed, yes_no},
};

/// Formats a twisted-Edwards curve compactly.
fn format_twisted_edwards_curve<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
where
    F::Elem: VisualizableField,
{
    let ax2 = if F::eq(curve.a(), &F::one()) {
        "x^2".to_string()
    } else {
        format!(
            "{}x^2",
            parenthesize_if_needed(&format_elem::<F>(curve.a()))
        )
    };
    let dx2y2 = if F::eq(curve.d(), &F::one()) {
        "x^2y^2".to_string()
    } else {
        format!(
            "{}x^2y^2",
            parenthesize_if_needed(&format_elem::<F>(curve.d()))
        )
    };

    format!("{} + y^2 = 1 + {}", ax2, dx2y2,)
}

/// Describes a twisted-Edwards curve in its native `a,d` presentation together
/// with the classical invariants derived from it.
fn describe_twisted_edwards_curve<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
where
    F::Elem: VisualizableField,
{
    [
        "Twisted-Edwards curve".to_string(),
        format!("equation: {}", format_twisted_edwards_curve(curve)),
        format!("characteristic: {}", F::characteristic()),
        format!("a: {}", format_elem::<F>(curve.a())),
        format!("d: {}", format_elem::<F>(curve.d())),
        "identity: (0, 1)".to_string(),
        "negation: -(x, y) = (-x, y)".to_string(),
        format!("discriminant: {}", format_elem::<F>(&curve.discriminant())),
        format!("c4: {}", format_elem::<F>(&curve.c4())),
        format!("c6: {}", format_elem::<F>(&curve.c6())),
        format!("j-invariant: {}", format_elem::<F>(&curve.j_invariant())),
    ]
    .join("\n")
}

/// Describes the canonical whole-curve bridge from the twisted-Edwards model
/// to its Montgomery companion.
fn describe_twisted_edwards_montgomery_companion<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let montgomery = curve.as_montgomery();

    [
        "Twisted-Edwards to Montgomery companion".to_string(),
        format!("source curve: {}", format_twisted_edwards_curve(curve)),
        format!("target curve: {}", format_montgomery_curve(&montgomery)),
        "curve-level formulas: A = 2(a + d)/(a - d), B = 4/(a - d)".to_string(),
        format!(
            "invariants preserved: c4={}, c6={}, discriminant={}, j={}",
            yes_no(F::eq(&curve.c4(), &montgomery.c4())),
            yes_no(F::eq(&curve.c6(), &montgomery.c6())),
            yes_no(F::eq(&curve.discriminant(), &montgomery.discriminant())),
            yes_no(F::eq(&curve.j_invariant(), &montgomery.j_invariant())),
        ),
        "point transport: total from affine twisted-Edwards points to Montgomery, but only partially defined in the reverse affine direction".to_string(),
    ]
    .join("\n")
}

/// Describes the current point transport between the affine twisted-Edwards
/// and Montgomery charts.
fn describe_twisted_edwards_birational_transport<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let montgomery = curve.as_montgomery();

    [
        "Birational point transport on the affine open".to_string(),
        format!(
            "twisted-Edwards chart: {}",
            format_twisted_edwards_curve(curve)
        ),
        format!("Montgomery chart: {}", format_montgomery_curve(&montgomery)),
        "current status: total in the direction Edwards -> Montgomery; still partial in the reverse affine direction".to_string(),
        "forward formulas: u = (1 + y)/(1 - y), v = (1 + y)/(x(1 - y))".to_string(),
        "inverse formulas: x = u/v, y = (u - 1)/(u + 1)".to_string(),
        "twisted-Edwards extension at x = 0: (0, 1) -> O and (0, -1) -> (0, 0)".to_string(),
        "reverse Montgomery exceptional locus: O, y = 0, or x = -1".to_string(),
        "status: this is not yet a total affine point correspondence in both directions".to_string(),
    ]
    .join("\n")
}

impl<F: Field> Visualizable for TwistedEdwardsCurve<F>
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    fn format_compact(&self) -> String {
        format_twisted_edwards_curve(self)
    }

    fn describe(&self) -> String {
        describe_twisted_edwards_curve(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::TwistedEdwardsCurve;
    use crate::fields::traits::Field;
    use crate::visualization::Visualizable;

    type F5 = crate::fields::Fp5;

    #[test]
    fn compact_formatter_shows_the_twisted_edwards_equation() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");

        assert_eq!(
            format_twisted_edwards_curve(&curve),
            "x^2 + y^2 = 1 + 2x^2y^2"
        );
    }

    #[test]
    fn curve_description_mentions_the_finite_identity_and_invariants() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");
        let description = describe_twisted_edwards_curve(&curve);

        assert!(description.contains("identity: (0, 1)"));
        assert!(description.contains("negation: -(x, y) = (-x, y)"));
        assert!(description.contains("j-invariant"));
    }

    #[test]
    fn montgomery_companion_description_mentions_curve_level_formulas() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");
        let description = describe_twisted_edwards_montgomery_companion(&curve);

        assert!(description.contains("A = 2(a + d)/(a - d), B = 4/(a - d)"));
        assert!(description.contains("total from affine twisted-Edwards points to Montgomery"));
    }

    #[test]
    fn birational_transport_description_mentions_the_current_extension_and_remaining_exceptions() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");
        let description = describe_twisted_edwards_birational_transport(&curve);

        assert!(description.contains("total in the direction Edwards -> Montgomery"));
        assert!(description.contains("(0, 1) -> O and (0, -1) -> (0, 0)"));
        assert!(description.contains("reverse Montgomery exceptional locus: O, y = 0, or x = -1"));
    }

    #[test]
    fn visualizable_impl_delegates_to_the_curve_helpers() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");

        assert_eq!(curve.format_compact(), "x^2 + y^2 = 1 + 2x^2y^2");
        assert!(curve.describe().contains("Twisted-Edwards curve"));
    }
}
