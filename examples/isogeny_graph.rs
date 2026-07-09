use crypto_bigint::{U64, const_prime_monty_params};
use num_bigint::BigUint;

use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::{
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    },
};
use elliptic_algorithms_lab::fields::{Fp, Fp7, Fp17};
use elliptic_algorithms_lab::isogenies::graphs::{
    IsogenyGraph, IsogenyGraphBuilder, IsogenyGraphNodeId,
    endomorphisms::refinement::CandidateRefinementStrategy,
};
use elliptic_algorithms_lab::visualization::Visualizable;

const_prime_monty_params!(Fp2749Params, U64, "0000000000000abd", 6);

type F = Fp17;
type F2749 = Fp<Fp2749Params, { U64::LIMBS }>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(1), F::from_i64(0))?;

    let graph = IsogenyGraphBuilder::new(curve.clone(), 2)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    println!("Educational ℓ-isogeny graph explorer");
    println!("======================================================");
    println!();
    println!("{}", graph.describe());
    println!();

    println!("Stored representative curves:");
    for (node_id, curve) in graph.node_representatives() {
        println!("  v{}: {curve}", node_id.0);
    }
    println!();

    let report = graph.verify_locally()?;
    println!("{}", report.describe());
    println!();

    let layers = graph.infer_volcano_like_layers(IsogenyGraphNodeId(0));
    println!("{}", layers.describe());
    println!();

    let endomorphism_report = graph.endomorphism_report_at(&BigUint::from(2u8))?;
    println!("{}", endomorphism_report.describe());
    println!();

    let refinement_report = endomorphism_report
        .refine_candidates_to_fixed_point(CandidateRefinementStrategy::default())?;
    println!("{}", refinement_report.describe());
    println!();

    let horizontal_curve = ShortWeierstrassCurve::<Fp7>::new(Fp7::from_i64(2), Fp7::from_i64(3))?;
    let horizontal_graph = IsogenyGraphBuilder::new(horizontal_curve, 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), BigUint::from(1u8))?;
    let prime_behavior = order.prime_behavior(&BigUint::from(3u8))?;
    let split_ideal = PrimeNormIdeal::split(order, BigUint::from(3u8), BigUint::from(1u8))?;
    let horizontal_ideal_reports =
        horizontal_graph.horizontal_ideal_reports(split_ideal.clone())?;
    let crater_walk_report =
        horizontal_graph.crater_walk_report(split_ideal, IsogenyGraphNodeId(0))?;

    println!("Horizontal ideal witnesses");
    println!("--------------------------");
    println!("Ideal side: in O_{{-23}}, 𝔭 = (3, ω - 1) has norm 3.");
    println!("Graph side: a 3-isogeny crater report supplies the horizontal-edge evidence.");
    println!();
    println!("{}", prime_behavior.describe());
    println!();
    println!(
        "Horizontal ideal compatibility reports: {}",
        horizontal_ideal_reports.len()
    );
    for report in &horizontal_ideal_reports {
        println!("{}", report.describe());
    }
    println!();
    println!("{}", crater_walk_report.describe());
    println!();

    let recovery_primes = primes(&[2, 3, 5]);
    let floor_graph = root_graph2749(curve2749(666, 215)?)?;
    let floor_partial =
        floor_graph.recover_endomorphism_ring_at(IsogenyGraphNodeId(0), &recovery_primes[..1])?;
    let floor_recovery =
        floor_graph.recover_endomorphism_ring_at(IsogenyGraphNodeId(0), &recovery_primes)?;

    println!("Multi-prime endomorphism-ring recovery example");
    println!("----------------------------------------------");
    println!("Curve over F_2749: y^2 = x^3 + 666x + 215");
    println!("For this curve, t = 14 and Δ_π = -10800 = 60^2 · (-3).");
    println!("So the Frobenius conductor is v = 60 = 2^2 · 3 · 5.");
    println!();
    for local in floor_recovery.local_reports() {
        println!("{}", local.describe());
        println!();
    }
    println!("Partial assembly from only the 2-volcano:");
    println!("{}", floor_partial.describe());
    println!();
    println!("Complete assembly from the 2-, 3-, and 5-volcanoes:");
    println!("{}", floor_recovery.describe());
    println!();

    let above_floor_graph = root_graph2749(curve2749(411, 1268)?)?;
    let above_floor_recovery =
        above_floor_graph.recover_endomorphism_ring_at(IsogenyGraphNodeId(0), &recovery_primes)?;
    let three_local = above_floor_recovery
        .local_reports()
        .iter()
        .find(|report| report.prime() == &BigUint::from(3u8))
        .expect("the recovery request included ℓ = 3");

    println!("A curve above the floor of the 3-volcano");
    println!("----------------------------------------");
    println!("Curve over F_2749: y^2 = x^3 + 411x + 1268");
    println!("This curve has the same Frobenius discriminant Δ_π = -10800.");
    println!("The 3-volcano probe certifies δ = 1, so the 3-part drops from O_60.");
    println!();
    println!("Local probe that sees the vertex is not on the floor:");
    println!("{}", three_local.describe());
    println!();
    println!("Complete assembly for this non-floor vertex:");
    println!("{}", above_floor_recovery.describe());

    Ok(())
}

fn curve2749(a: i64, b: i64) -> Result<ShortWeierstrassCurve<F2749>, Box<dyn std::error::Error>> {
    Ok(ShortWeierstrassCurve::<F2749>::new(
        F2749::from_i64(a),
        F2749::from_i64(b),
    )?)
}

fn root_graph2749(
    curve: ShortWeierstrassCurve<F2749>,
) -> Result<IsogenyGraph<ShortWeierstrassCurve<F2749>>, Box<dyn std::error::Error>> {
    Ok(IsogenyGraphBuilder::new(curve, 2)
        .max_depth(0)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?)
}

fn primes(values: &[usize]) -> Vec<BigUint> {
    values.iter().copied().map(BigUint::from).collect()
}
