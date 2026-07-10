use num_bigint::BigUint;

use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::{
        binary_quadratic_forms::QuadraticClassGroup,
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    },
};
use elliptic_algorithms_lab::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};
use elliptic_algorithms_lab::visualization::Visualizable;

type F7 = elliptic_algorithms_lab::fields::Fp7;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))?;
    let graph = IsogenyGraphBuilder::new(curve.clone(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), BigUint::from(1u8))?;
    let ideal = PrimeNormIdeal::split(order, BigUint::from(3u8), BigUint::from(1u8))?;
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))?;
    let start = IsogenyGraphNodeId(0);

    let crater = graph.volcano_crater_report(ideal.norm())?;
    let labeled_walk = graph.labeled_crater_walk_report(&class_group, ideal, start)?;

    println!("Crater walk labeled by an ideal/form class");
    println!("==========================================");
    println!();
    println!("Setup");
    println!("-----");
    println!("Curve: {} over F_7", curve.format_compact());
    println!("Local isogeny degree: ℓ = 3");
    println!("Quadratic order: discriminant D = -23");
    println!("Prime ideal: 𝔭 = (3, ω - 1)");
    println!();

    println!("{}", crater.describe());
    println!();

    println!("{}", labeled_walk.describe());
    println!();
    println!("What this certifies");
    println!("-------------------");
    println!("The ideal, the reduced form class, and the local crater prime are compatible.");
    println!("The recorded walk follows certified horizontal crater edges in graph order.");

    Ok(())
}
