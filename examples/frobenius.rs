use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    frobenius::{
        characteristic_equation::FrobeniusCharacteristicEquationCurveModel,
        extension_counts::compare_extension_count_with_enumeration,
    },
    short_weierstrass::isomorphisms::ShortWeierstrassQuadraticTwist,
    traits::{AffineCurveModel, EnumerableCurveModel, FrobeniusTraceCurveModel},
};
use elliptic_algorithms_lab::fields::traits::*;
use elliptic_algorithms_lab::fields::traits::{EnumerableFiniteField, SqrtField};
use elliptic_algorithms_lab::isogenies::{
    frobenius_relation::FrobeniusComparableIsogeny,
    scalar_multiplication::ScalarMultiplicationIsogeny,
};
use elliptic_algorithms_lab::visualization::{
    Visualizable, describe_frobenius_characteristic_equation_check,
    describe_frobenius_extension_enumeration_comparison_report, describe_hasse_bound_report,
    describe_isogeny_frobenius_relation, describe_quadratic_twist_frobenius_relation, format_curve,
    format_frobenius_trace, format_point_compact,
};

type F17 = elliptic_algorithms_lab::fields::Fp17;
type F19 = elliptic_algorithms_lab::fields::Fp19;
type F41 = elliptic_algorithms_lab::fields::Fp41;
type F43 = elliptic_algorithms_lab::fields::Fp43;

elliptic_algorithms_lab::fields::extension_field::define_fp_quadratic_extension!(
    spec: ExampleF17Sqrt3Spec,
    field: ExampleF17Sqrt3,
    base: F17,
    non_residue: 3,
    name: "ExampleF17(sqrt(3))",
);

fn heading(title: &str) {
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

fn first_nonsquare<F>() -> F::Elem
where
    F: EnumerableFiniteField + SqrtField,
{
    F::elements()
        .into_iter()
        .find(|value| !F::is_zero(value) && !F::has_square_root(value))
        .expect("small odd prime fields should contain non-squares")
}

fn lift_f17_curve_to_quadratic_extension(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<ExampleF17Sqrt3> {
    ShortWeierstrassCurve::<ExampleF17Sqrt3>::new(
        ExampleF17Sqrt3::from_base(*curve.a()),
        ExampleF17Sqrt3::from_base(*curve.b()),
    )
    .expect("lifting a smooth F17 curve to F17^2 should preserve smoothness")
}

fn first_non_fixed_extension_point(
    curve: &ShortWeierstrassCurve<ExampleF17Sqrt3>,
) -> Option<AffinePoint<ExampleF17Sqrt3>> {
    curve.points().into_iter().find(|point| {
        curve
            .absolute_frobenius_orbit(point, 1)
            .map(|orbit| orbit.period() > 1)
            .unwrap_or(false)
    })
}

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Frobenius milestone tour");
    println!("======================================================");
    println!();

    heading("1. Prime-field package");
    let prime_curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one())?;
    let prime_trace = prime_curve.frobenius_trace()?;
    let prime_polynomial = prime_trace.characteristic_polynomial();
    let prime_zeta = prime_trace.local_zeta_function();
    let hasse_report = prime_curve.verify_hasse_bound()?;
    let curve_type = prime_trace.curve_type();

    println!("curve: {}", format_curve(&prime_curve));
    println!("trace: {}", format_frobenius_trace(&prime_trace));
    println!("χ_π(T): {}", prime_polynomial);
    println!("Z(E/F_q, T): {}", prime_zeta);
    println!("curve type: {}", curve_type.format_compact());
    println!("{}", indent(&describe_hasse_bound_report(&hasse_report), 2));
    println!();

    heading("2. Characteristic equation at one point");
    let rational_point = prime_curve.point(F43::zero(), F43::one())?;
    let check = prime_curve
        .verify_frobenius_characteristic_equation_at_point(&rational_point, &prime_polynomial)?;

    println!("chosen point: {}", format_point_compact(&rational_point));
    println!(
        "{}",
        indent(&describe_frobenius_characteristic_equation_check(&check), 2)
    );
    println!();

    heading("3. Extension-field viewpoint");
    let base_curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))?;
    let extension_curve = lift_f17_curve_to_quadratic_extension(&base_curve);
    let base_trace = base_curve.frobenius_trace()?;
    let extension_comparison =
        compare_extension_count_with_enumeration(&extension_curve, &base_trace)?;

    println!("base curve over F_17: {}", format_curve(&base_curve));
    println!(
        "{}",
        indent(
            &describe_frobenius_extension_enumeration_comparison_report(&extension_comparison),
            2,
        )
    );

    if let Some(point) = first_non_fixed_extension_point(&extension_curve) {
        let orbit = extension_curve.absolute_frobenius_orbit(&point, 1)?;
        println!("sample point visible only over the quadratic extension:");
        println!("  start = {}", point.format_compact());
        println!("  orbit = {}", orbit.format_compact());
    } else {
        println!("all enumerated points happened to be fixed by π_17 in this sample.");
    }
    println!();

    heading("4. Twist and isogeny relations");
    let twist_curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))?;
    let twist_package = ShortWeierstrassQuadraticTwist::new(twist_curve, first_nonsquare::<F19>())?;
    let twist_relation = twist_package.frobenius_relation()?;

    let scalar_isogeny = ScalarMultiplicationIsogeny::new(f41_curve(), 2)?;
    let isogeny_relation = scalar_isogeny.frobenius_relation_report()?;

    println!(
        "{}",
        indent(
            &describe_quadratic_twist_frobenius_relation(&twist_relation),
            2,
        )
    );
    println!();
    println!(
        "{}",
        indent(&describe_isogeny_frobenius_relation(&isogeny_relation), 2)
    );

    println!();
    println!("Conclusions");
    println!("-----------");
    println!(
        "  1. The trace package is the central hub: χ_π(T), Z(E/F_q, T), Hasse, and the curve type all come from it."
    );
    println!(
        "  2. Over the represented base field, π_q behaves trivially on rational points, but π_p over an extension already detects larger fields of definition."
    );
    println!(
        "  3. Point-count relations for twists and isogenies become easy to inspect once everything is phrased in terms of Frobenius traces."
    );

    Ok(())
}
