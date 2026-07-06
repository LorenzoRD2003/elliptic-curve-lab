use crypto_bigint::{U64, const_prime_monty_params};
use num_bigint::BigUint;

use elliptic_algorithms_lab::elliptic_curves::ShortWeierstrassCurve;
use elliptic_algorithms_lab::fields::{Fp, Fp17};
use elliptic_algorithms_lab::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId, endomorphisms::refinement::CandidateRefinementStrategy,
};
use elliptic_algorithms_lab::visualization::{
    EndomorphismRingRecoveryAssembly as Assembly, explain_graph_candidate_refinement_report,
    explain_graph_endomorphism_report, explain_isogeny_graph,
    explain_short_weierstrass_root_endomorphism_ring_level_recovery_walkthrough as explain_ring_recovery,
    explain_volcano_like_layers,
};

const_prime_monty_params!(Fp2749Params, U64, "0000000000000abd", 6);

type F = Fp17;
type F2749 = Fp<Fp2749Params, { U64::LIMBS }>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(1), F::from_i64(0))?;

    let graph = IsogenyGraphBuilder::new(curve, 2)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    println!("Educational ℓ-isogeny graph explorer");
    println!("======================================================");
    println!();
    println!("{}", explain_isogeny_graph(&graph));
    println!();

    let report = graph.verify_locally()?;
    println!("Local verification report:");
    println!("{report:#?}");
    println!();

    let layers = graph.infer_volcano_like_layers(IsogenyGraphNodeId(0));
    println!("{}", explain_volcano_like_layers(&graph, &layers));
    println!();

    let endomorphism_report = graph.endomorphism_report_at(&BigUint::from(2u8))?;
    println!(
        "{}",
        explain_graph_endomorphism_report(&endomorphism_report)
    );
    println!();

    let refinement_report = endomorphism_report
        .refine_candidates_to_fixed_point(CandidateRefinementStrategy::default())?;
    println!(
        "{}",
        explain_graph_candidate_refinement_report(&refinement_report)
    );
    println!();

    println!(
        "{}",
        explain_ring_recovery(
            "Multi-prime endomorphism-ring recovery example",
            &[
                "Curve over F_2749: y^2 = x^3 + 666x + 215",
                "For this curve, t = 14 and Δ_π = -10800 = 60^2 · (-3).",
                "So the Frobenius conductor is v = 60 = 2^2 · 3 · 5.",
            ],
            curve2749(666, 215)?,
            &[2, 3, 5],
            &[("", 2), ("", 3), ("", 5)],
            &[
                Assembly::first_primes("Partial assembly from only the 2-volcano:", 1),
                Assembly::complete("Complete assembly from the 2-, 3-, and 5-volcanoes:"),
            ],
        )?
    );
    println!();

    println!(
        "{}",
        explain_ring_recovery(
            "A curve above the floor of the 3-volcano",
            &[
                "Curve over F_2749: y^2 = x^3 + 411x + 1268",
                "This curve has the same Frobenius discriminant Δ_π = -10800.",
                "The 3-volcano probe certifies δ = 1, so the 3-part drops from O_60.",
            ],
            curve2749(411, 1268)?,
            &[2, 3, 5],
            &[("Local probe that sees the vertex is not on the floor:", 3)],
            &[Assembly::complete(
                "Complete assembly for this non-floor vertex:",
            )],
        )?
    );

    Ok(())
}

fn curve2749(a: i64, b: i64) -> Result<ShortWeierstrassCurve<F2749>, Box<dyn std::error::Error>> {
    Ok(ShortWeierstrassCurve::<F2749>::new(
        F2749::from_i64(a),
        F2749::from_i64(b),
    )?)
}
