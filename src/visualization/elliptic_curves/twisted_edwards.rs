use core::fmt;

use crate::elliptic_curves::TwistedEdwardsCurve;
use crate::fields::traits::Field;
use crate::visualization::{
    elliptic_curves::montgomery::format_montgomery_curve, fields::traits::VisualizableField,
    traits::Visualizable,
};

fn format_elem<F: Field>(value: &F::Elem) -> String
where
    F::Elem: VisualizableField,
{
    value.format_elem()
}

fn parenthesize_if_needed(text: &str) -> String {
    if text.contains(' ') || text.starts_with('-') || text.contains('/') {
        format!("({text})")
    } else {
        text.to_string()
    }
}

/// Formats a twisted-Edwards curve compactly.
pub fn format_twisted_edwards_curve<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
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
pub fn describe_twisted_edwards_curve<F: Field>(curve: &TwistedEdwardsCurve<F>) -> String
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
pub fn describe_twisted_edwards_montgomery_companion<F: Field>(
    curve: &TwistedEdwardsCurve<F>,
) -> String
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
            if F::eq(&curve.c4(), &montgomery.c4()) {
                "yes"
            } else {
                "no"
            },
            if F::eq(&curve.c6(), &montgomery.c6()) {
                "yes"
            } else {
                "no"
            },
            if F::eq(&curve.discriminant(), &montgomery.discriminant()) {
                "yes"
            } else {
                "no"
            },
            if F::eq(&curve.j_invariant(), &montgomery.j_invariant()) {
                "yes"
            } else {
                "no"
            },
        ),
        "point transport: only partial at this stage, through the birational affine-open formulas"
            .to_string(),
    ]
    .join("\n")
}

/// Describes the current partial point transport between the affine
/// twisted-Edwards and Montgomery charts.
pub fn describe_twisted_edwards_birational_transport<F: Field>(
    curve: &TwistedEdwardsCurve<F>,
) -> String
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
        "forward formulas: u = (1 + y)/(1 - y), v = (1 + y)/(x(1 - y))".to_string(),
        "inverse formulas: x = u/v, y = (u - 1)/(u + 1)".to_string(),
        "twisted-Edwards exceptional locus: x = 0 or y = 1".to_string(),
        "Montgomery exceptional locus: O, y = 0, or x = -1".to_string(),
        "status: this is a partial birational transport, not a total affine point correspondence"
            .to_string(),
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
    use crate::elliptic_curves::TwistedEdwardsCurve;
    use crate::fields::{Fp, traits::Field};
    use crate::visualization::{
        elliptic_curves::twisted_edwards::{
            describe_twisted_edwards_birational_transport, describe_twisted_edwards_curve,
            describe_twisted_edwards_montgomery_companion, format_twisted_edwards_curve,
        },
        traits::Visualizable,
    };

    type F5 = Fp<5>;

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
        assert!(description.contains("point transport: only partial"));
    }

    #[test]
    fn birational_transport_description_mentions_the_exceptional_loci() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");
        let description = describe_twisted_edwards_birational_transport(&curve);

        assert!(description.contains("twisted-Edwards exceptional locus: x = 0 or y = 1"));
        assert!(description.contains("Montgomery exceptional locus: O, y = 0, or x = -1"));
        assert!(description.contains("partial birational transport"));
    }

    #[test]
    fn visualizable_impl_delegates_to_the_curve_helpers() {
        let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
            .expect("sample twisted-Edwards curve should be non-singular");

        assert_eq!(curve.format_compact(), "x^2 + y^2 = 1 + 2x^2y^2");
        assert!(curve.describe().contains("Twisted-Edwards curve"));
    }
}
