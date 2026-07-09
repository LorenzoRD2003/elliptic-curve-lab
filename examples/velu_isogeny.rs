use elliptic_algorithms_lab::elliptic_curves::short_weierstrass::isogenies::VeluIsogeny;
use elliptic_algorithms_lab::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, traits::AffineCurveModel,
};
use elliptic_algorithms_lab::isogenies::traits::Isogeny;
use elliptic_algorithms_lab::visualization::Visualizable;

type F = elliptic_algorithms_lab::fields::Fp101;

fn main() -> Result<(), CurveError> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let generator = curve.point(F::from_i64(35), F::from_i64(15))?;
    let point = curve.point(F::from_i64(1), F::from_i64(39))?;

    let isogeny =
        VeluIsogeny::from_generator(curve.clone(), generator.clone()).expect("Vélu isogeny");
    let image = isogeny
        .evaluate(&point)
        .expect("the chosen point should map into the codomain");

    println!("Vélu isogeny over a small prime field");
    println!("=======================================================");
    println!();
    println!("{}", isogeny.describe());
    println!();
    println!("sample evaluation:");
    println!("  P = {}", generator.format_compact());
    println!("  Q = {}", point.format_compact());
    println!("  phi(Q) = {}", image.format_compact());
    println!();
    println!("codomain:");
    println!("{}", isogeny.codomain().describe());

    Ok(())
}
