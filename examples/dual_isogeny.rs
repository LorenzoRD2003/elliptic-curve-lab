use elliptic_algorithms_lab::{
    AffineCurveModel, CurveError, Field, Fp, Isogeny, ShortWeierstrassCurve, VeluIsogeny,
    describe_dual_isogeny, format_curve, format_point_compact, summarize_dual_verification,
};

type F = Fp<29>;

fn indent_block(block: &str) -> String {
    block
        .lines()
        .map(|line| format!("  {line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn main() -> Result<(), CurveError> {
    let domain = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(2))?;
    let generator = domain.point(F::from_i64(10), F::from_i64(23))?;
    let phi = VeluIsogeny::from_generator(domain.clone(), generator.clone()).expect("Vélu isogeny");
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found by exhaustive search");

    println!("Fifth milestone: dual isogenies and composition");
    println!("===============================================");
    println!();
    println!("domain E:");
    println!("  {}", format_curve(&domain));
    println!();
    println!("isogeny phi:");
    println!("  degree: {}", phi.degree());
    println!("  kernel: <P>");
    println!("  P = {}", format_point_compact(&generator));
    println!();
    println!("codomain E':");
    println!("  {}", format_curve(phi.codomain()));
    println!();
    println!("dual phi_hat:");
    println!("{}", indent_block(&describe_dual_isogeny(&dual)));
    println!();
    println!("checks:");
    println!(
        "{}",
        indent_block(
            &summarize_dual_verification(&phi, &dual)
                .expect("dual verification summary should build"),
        )
    );

    Ok(())
}
