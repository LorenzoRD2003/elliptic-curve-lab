use elliptic_algorithms_lab::elliptic_curves::ShortWeierstrassQuadraticTwist;
use elliptic_algorithms_lab::{
    Field, Fp, ShortWeierstrassCurve, describe_isomorphism, explain_quadratic_twist,
    explain_short_weierstrass_scaling, format_extension_field, format_extension_field_element,
    summarize_curve_comparison,
};

type F7 = Fp<7>;
type F13 = Fp<13>;
type F19 = Fp<19>;

elliptic_algorithms_lab::fields::define_fp_quadratic_extension!(
    spec: F19Sqrt2Spec,
    field: F19Sqrt2,
    base: F19,
    non_residue: 2,
    name: "F19(sqrt(2))",
);

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Isomorphisms, twists, and short-Weierstrass comparison");
    println!("====================================================================");
    println!();

    scaling_example()?;
    quadratic_twist_example()?;
    automorphism_example()?;

    Ok(())
}

fn scaling_example() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))?;
    let scaled = curve.scaled_by(F7::from_i64(3))?;
    let isomorphism = curve
        .find_isomorphism_to(&scaled)
        .expect("a scaled curve should be found by exhaustive search over F7");

    println!("1. Same curve up to a change of coordinates");
    println!("-------------------------------------------");
    println!(
        "{}",
        explain_short_weierstrass_scaling(&curve, &F7::from_i64(3))?
    );
    println!();
    println!("{}", describe_isomorphism(&isomorphism));
    println!();
    println!("{}", summarize_curve_comparison(&curve, &scaled));
    println!();

    Ok(())
}

fn quadratic_twist_example() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))?;
    let trivial_factor = F19::from_i64(4);
    let quadratic_factor = F19::from_i64(2);
    let trivial_twist = curve.quadratic_twist(trivial_factor)?;
    let quadratic_twist = curve.quadratic_twist(quadratic_factor)?;
    let package = ShortWeierstrassQuadraticTwist::new(curve.clone(), quadratic_factor)?;
    let extension_isomorphism = package.isomorphism_over_quadratic_extension::<F19Sqrt2Spec>()?;

    println!("2. Same j-invariant versus base-field isomorphism");
    println!("-------------------------------------------------");
    println!("Twist by a square:");
    println!("{}", explain_quadratic_twist(&curve, &trivial_factor)?);
    println!();
    println!("{}", summarize_curve_comparison(&curve, &trivial_twist));
    println!();
    println!("Twist by a non-square:");
    println!("{}", explain_quadratic_twist(&curve, &quadratic_factor)?);
    println!();
    println!("{}", summarize_curve_comparison(&curve, &quadratic_twist));
    println!();
    println!("Over the quadratic extension:");
    println!(
        "  ambient field: {}",
        format_extension_field::<F19Sqrt2Spec>()
    );
    println!(
        "  scaling witness u: {}",
        format_extension_field_element::<F19Sqrt2Spec>(extension_isomorphism.scaling_factor())
    );
    println!("  check: u^2 = d inside the extension, so the twist becomes isomorphic there");
    println!();

    Ok(())
}

fn automorphism_example() -> Result<(), Box<dyn std::error::Error>> {
    let generic = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))?;
    let j_1728 = ShortWeierstrassCurve::<F13>::new(F13::from_i64(1), F13::zero())?;
    let j_zero = ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::one())?;

    println!("3. Automorphisms and special j-loci");
    println!("-----------------------------------");
    println!(
        "Generic curve over F7: {} automorphisms",
        generic.automorphisms().len()
    );
    println!("  expected: only +/-1");
    for iso in generic.automorphisms() {
        println!("  u = {}", iso.scaling_factor());
    }
    println!();
    println!(
        "j = 1728 curve over F13: {} automorphisms",
        j_1728.automorphisms().len()
    );
    println!("  this is the b = 0 locus, so extra automorphisms can appear when i is available");
    for iso in j_1728.automorphisms() {
        println!("  u = {}", iso.scaling_factor());
    }
    println!();
    println!(
        "j = 0 curve over F13: {} automorphisms",
        j_zero.automorphisms().len()
    );
    println!(
        "  this is the a = 0 locus, so extra automorphisms can appear when non-trivial cube roots are available"
    );
    for iso in j_zero.automorphisms() {
        println!("  u = {}", iso.scaling_factor());
    }
    println!();

    Ok(())
}
