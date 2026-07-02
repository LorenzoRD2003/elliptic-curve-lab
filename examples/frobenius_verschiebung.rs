use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::function_fields::ShortWeierstrassFunction,
    short_weierstrass::isogenies::frobenius::FrobeniusVerschiebungFactorizationReport,
    short_weierstrass::isogenies::function_field_maps::ShortWeierstrassFunctionFieldMap,
    traits::EnumerableCurveModel,
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::fields::{
    FieldError,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
    rational_function_field::RationalFunction,
};
use elliptic_algorithms_lab::isogenies::scalar_multiplication::ScalarMultiplicationIsogeny;
use elliptic_algorithms_lab::isogenies::traits::Isogeny;
use elliptic_algorithms_lab::visualization::*;
use elliptic_algorithms_lab::visualization::{
    Visualizable, describe_differential_pullback_report,
    describe_frobenius_verschiebung_factorization_report,
    explain_frobenius_verschiebung_factorization_report, fields::describe_extension_field,
    format_curve, format_point_compact, format_short_weierstrass_function_field_map,
};
use num_traits::ToPrimitive;

type F5 = elliptic_algorithms_lab::fields::Fp5;

#[derive(Clone, Copy)]
struct F25ExampleSpec;

impl ExtensionFieldSpec for F25ExampleSpec {
    type Base = F5;

    const NAME: &'static str = "F_25 example field";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![
            Self::Base::from_i64(-2),
            Self::Base::zero(),
            Self::Base::one(),
        ])
        .expect("x^2 - 2 should be a valid structural modulus over F5")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

type F25Example = ExtensionField<F25ExampleSpec>;
type SamplePointStoryRow = (
    AffinePoint<F25Example>,
    AffinePoint<F25Example>,
    AffinePoint<F25Example>,
    AffinePoint<F25Example>,
);

fn heading(title: &str) {
    println!("{title}");
    println!("{}", "=".repeat(title.len()));
}

fn section(title: &str) {
    println!("{title}");
    println!("{}", "-".repeat(title.len()));
}

fn indent(block: &str, spaces: usize) -> String {
    let prefix = " ".repeat(spaces);
    block
        .lines()
        .map(|line| format!("{prefix}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn alpha() -> <F25Example as Field>::Elem {
    F25Example::element(vec![F5::zero(), F5::one()])
}

fn curve() -> ShortWeierstrassCurve<F25Example> {
    ShortWeierstrassCurve::new(alpha(), F25Example::one())
        .expect("valid nontrivial curve over F_25")
}

fn evaluate_rational_function_at_x(
    function: &RationalFunction<F25Example>,
    x: &<F25Example as Field>::Elem,
) -> Option<<F25Example as Field>::Elem> {
    let numerator = function.numerator().evaluate(x).ok()?;
    let denominator = function.denominator().evaluate(x).ok()?;

    if F25Example::is_zero(&denominator) {
        None
    } else {
        F25Example::div(&numerator, &denominator).ok()
    }
}

fn evaluate_function_at_point(
    function: &ShortWeierstrassFunction<F25Example>,
    point: &AffinePoint<F25Example>,
) -> Option<<F25Example as Field>::Elem> {
    match point {
        AffinePoint::Infinity => None,
        AffinePoint::Finite { x, y } => {
            let a_value = evaluate_rational_function_at_x(function.a_part(), x)?;
            let b_value = evaluate_rational_function_at_x(function.b_part(), x)?;
            Some(F25Example::add(&a_value, &F25Example::mul(y, &b_value)))
        }
    }
}

fn evaluate_map_at_point(
    map: &ShortWeierstrassFunctionFieldMap<F25Example>,
    point: &AffinePoint<F25Example>,
) -> Option<AffinePoint<F25Example>> {
    match point {
        AffinePoint::Infinity => Some(AffinePoint::Infinity),
        AffinePoint::Finite { .. } => {
            let x = evaluate_function_at_point(map.x_pullback(), point)?;
            let y = evaluate_function_at_point(map.y_pullback(), point)?;
            Some(AffinePoint::new(x, y))
        }
    }
}

fn sample_point_story(
    curve: &ShortWeierstrassCurve<F25Example>,
    report: &FrobeniusVerschiebungFactorizationReport<F25Example>,
    multiplication_by_p: &ScalarMultiplicationIsogeny<ShortWeierstrassCurve<F25Example>>,
) -> Vec<SamplePointStoryRow> {
    curve
        .points()
        .into_iter()
        .filter(|point| matches!(point, AffinePoint::Finite { .. }))
        .filter_map(|point| {
            let frobenius_image = report.frobenius().evaluate(&point).ok()?;
            let v_after_f = evaluate_map_at_point(
                report.verschiebung().as_function_field_map(),
                &frobenius_image,
            )?;
            let scalar_image = multiplication_by_p.evaluate(&point).ok()?;

            Some((point, frobenius_image, v_after_f, scalar_image))
        })
        .take(3)
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = curve();
    let multiplication_by_p = ScalarMultiplicationIsogeny::new(
        curve.clone(),
        F25Example::characteristic()
            .to_positive_biguint()
            .and_then(|p| p.to_u64())
            .expect("example characteristic should fit the scalar API"),
    )?;
    let report = multiplication_by_p.frobenius_verschiebung_factorization_report()?;
    let certified_p_pullback =
        multiplication_by_p.as_function_field_map_from_verschiebung(report.certificate())?;
    let point_stories = sample_point_story(&curve, &report, &multiplication_by_p);

    heading("Frobenius, Verschiebung, and [p] over F_(5^2)");
    println!();

    section("0. Ambient field and genuinely moved twist");
    println!(
        "{}",
        indent(&describe_extension_field::<F25ExampleSpec>(), 2)
    );
    println!();
    println!("chosen curve E: {}", format_curve(&curve));
    println!(
        "chosen coefficients: a = {}, b = {}",
        curve.a().format_elem(),
        curve.b().format_elem()
    );
    println!(
        "Frobenius moves the non-prime-field coefficient: a^5 = {}",
        F25Example::pow(
            curve.a(),
            &F25Example::characteristic()
                .to_positive_biguint()
                .expect("finite fields have positive characteristic"),
        )
        .format_elem()
    );
    println!(
        "codomain curve E^(5): {}",
        format_curve(report.frobenius().codomain())
    );
    println!(
        "is the Frobenius codomain literally the same displayed equation as E? {}",
        if report.frobenius().codomain() == &curve {
            "yes"
        } else {
            "no"
        }
    );
    println!(
        "scope note: this example studies the absolute Frobenius factorization [p] = V ∘ Frob_p; it is not a relative-Frobenius example."
    );
    println!();

    section("1. Compact summary");
    println!("{}", report.format_compact());
    println!();
    println!(
        "{}",
        indent(
            &describe_frobenius_verschiebung_factorization_report(&report),
            2,
        )
    );
    println!();

    section("2. The actual factorization story");
    println!(
        "{}",
        indent(
            &explain_frobenius_verschiebung_factorization_report(&report),
            2,
        )
    );
    println!();

    section("3. Direct [p]^* versus certified [p]^*");
    println!("curve E: {}", format_curve(&curve));
    println!(
        "absolute Frobenius codomain: {}",
        format_curve(report.frobenius().codomain())
    );
    println!(
        "direct [p]^*: {}",
        format_short_weierstrass_function_field_map(report.multiplication_by_p_pullback())
    );
    println!(
        "certified [p]^* from V and Frobenius: {}",
        format_short_weierstrass_function_field_map(&certified_p_pullback)
    );
    println!(
        "direct and certified pullbacks agree: {}",
        certified_p_pullback == *report.multiplication_by_p_pullback()
    );
    println!();

    section("4. Differential viewpoint");
    let frobenius_differential = report.frobenius_differential_report()?;
    let verschiebung_differential = report.verschiebung_differential_report()?;

    println!("Frobenius:");
    println!(
        "{}",
        indent(
            &describe_differential_pullback_report(&frobenius_differential),
            2,
        )
    );
    println!();
    println!("Verschiebung:");
    println!(
        "{}",
        indent(
            &describe_differential_pullback_report(&verschiebung_differential),
            2,
        )
    );
    println!();

    section("5. Sample points");
    for (index, (point, frobenius_image, v_after_f, scalar_image)) in
        point_stories.iter().enumerate()
    {
        println!("sample {}:", index + 1);
        println!("  P = {}", format_point_compact(point));
        println!("  Frob_p(P) = {}", format_point_compact(frobenius_image));
        println!("  V(Frob_p(P)) = {}", format_point_compact(v_after_f));
        println!("  [p]P = {}", format_point_compact(scalar_image));
        println!("  agreement: {}", v_after_f == scalar_image);
        println!();
    }

    section("6. Explicit verification");
    report.verify()?;
    report.certificate().verify_v_after_f_equals_p()?;
    report.certificate().verify_f_after_v_equals_p()?;
    println!("Both identities were rechecked successfully:");
    println!("  V ∘ Frob_p = [p]_E");
    println!("  Frob_p ∘ V = [p]_{{E^(p)}}");
    println!();

    section("Conclusions");
    println!(
        "1. Over F_(5^2), absolute Frobenius visibly moves the curve: E and E^(5) have different displayed equations because the coefficient a is not fixed by x ↦ x^5."
    );
    println!(
        "2. Even in that genuinely nontrivial situation, the direct generic-point pullback [5]^* still contains enough information to reconstruct Verschiebung."
    );
    println!(
        "3. The differential reports keep the expected arithmetic meaning: Frobenius is purely inseparable, while Verschiebung has non-zero differential multiplier and is therefore detected as separable."
    );
    println!(
        "4. The certified route and the direct route to [5]^* agree exactly, which is the main consistency check behind the current implementation."
    );

    Ok(())
}
