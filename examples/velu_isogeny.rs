use elliptic_algorithms_lab::{
    AffineCurveModel, CurveError, Field, Fp, Isogeny, ShortWeierstrassCurve, VeluIsogeny,
    describe_isogeny, explain_velu_codomain, explain_velu_evaluation, format_point_compact,
};

type F = Fp<101>;

fn main() -> Result<(), CurveError> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let generator = curve.point(F::from_i64(35), F::from_i64(15))?;
    let point = curve.point(F::from_i64(1), F::from_i64(39))?;

    let isogeny =
        VeluIsogeny::from_generator(curve.clone(), generator.clone()).expect("Vélu isogeny");
    let image = isogeny
        .evaluate(&point)
        .expect("the chosen point should map into the codomain");

    println!("Third milestone: Vélu isogeny over a small prime field");
    println!("=======================================================");
    println!();
    println!("{}", describe_isogeny(&isogeny));
    println!();
    println!("sample evaluation:");
    println!("  P = {}", format_point_compact(&generator));
    println!("  Q = {}", format_point_compact(&point));
    println!("  phi(Q) = {}", format_point_compact(&image));
    println!();
    println!("codomain derivation:");
    println!("{}", explain_velu_codomain(&isogeny));
    println!();
    println!("point evaluation derivation:");
    println!(
        "{}",
        explain_velu_evaluation(&isogeny, &point).expect("evaluation explanation")
    );

    Ok(())
}
