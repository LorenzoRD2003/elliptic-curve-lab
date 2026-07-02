use elliptic_algorithms_lab::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, traits::EnumerableCurveModel,
};
use elliptic_algorithms_lab::visualization::{format_curve, format_point};

type F = elliptic_algorithms_lab::fields::Fp101;

fn main() -> Result<(), CurveError> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let points = curve.points();
    let order = curve.order();

    println!("Curve order over a small prime field");
    println!("====================================================");
    println!();
    println!("curve: {}", format_curve(&curve));
    println!("#E(F_101) = {order}");
    println!("sample of enumerated points:");

    for point in points.iter().take(10) {
        println!("  {}", format_point(point));
    }

    if points.len() > 10 {
        println!("  ...");
    }

    Ok(())
}
