use elliptic_algorithms_lab::elliptic_curves::short_weierstrass::isogenies::{
    DualVeluIsogeny, VeluIsogeny,
};
use elliptic_algorithms_lab::elliptic_curves::{ShortWeierstrassCurve, traits::AffineCurveModel};
use elliptic_algorithms_lab::isogenies::traits::Isogeny;
use elliptic_algorithms_lab::visualization::Visualizable;

type F = elliptic_algorithms_lab::fields::Fp29;

fn indent_block(block: &str) -> String {
    block
        .lines()
        .map(|line| format!("  {line}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let domain = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(2))?;
    let generator = domain.point(F::from_i64(10), F::from_i64(23))?;
    let phi = VeluIsogeny::from_generator(domain.clone(), generator.clone()).expect("Vélu isogeny");
    let dual = phi
        .find_dual_exhaustively()
        .expect("dual should be found by exhaustive search");

    println!("Dual isogenies and composition");
    println!("===============================================");
    println!();
    println!("domain E:");
    println!("  {}", domain.format_compact());
    println!();
    println!("isogeny phi:");
    println!("  degree: {}", phi.degree());
    println!("  kernel: <P>");
    println!("  P = {}", generator.format_compact());
    println!();
    println!("codomain E':");
    println!("  {}", phi.codomain().format_compact());
    println!();
    println!("dual phi_hat:");
    println!("{}", indent_block(&dual.describe()));
    println!();
    println!("checks:");
    let report = DualVeluIsogeny::dual_report(&phi, &dual)?;
    println!("{}", indent_block(&report.describe()));

    Ok(())
}
