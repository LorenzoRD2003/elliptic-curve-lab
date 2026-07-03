use num_bigint::BigInt;
use num_rational::BigRational;

use elliptic_algorithms_lab::elliptic_curves::ShortWeierstrassCurve;
use elliptic_algorithms_lab::fields::Q;
use elliptic_algorithms_lab::visualization::describe_rational_torsion_report;

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn show_case(
    title: &str,
    curve: ShortWeierstrassCurve<Q>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("{title}");
    println!("{}", "-".repeat(title.len()));
    let report = curve.rational_torsion()?;
    println!("{}", describe_rational_torsion_report(&report));
    println!();
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Rational torsion for short-Weierstrass curves over Q");
    println!("====================================================");
    println!();
    println!(
        "Route: scale to an integral short-Weierstrass model, enumerate \
         Lutz-Nagell candidates, then verify exact Mazur point orders."
    );
    println!();

    show_case(
        "trivial torsion: y^2 = x^3 + x + 1",
        ShortWeierstrassCurve::<Q>::new(q(1, 1), q(1, 1))?,
    )?;

    show_case(
        "cyclic torsion: y^2 = x^3 + 1",
        ShortWeierstrassCurve::<Q>::new(q(0, 1), q(1, 1))?,
    )?;

    show_case(
        "product torsion after scaling: y^2 = x^3 - x/16",
        ShortWeierstrassCurve::<Q>::new(q(-1, 16), q(0, 1))?,
    )?;

    Ok(())
}
