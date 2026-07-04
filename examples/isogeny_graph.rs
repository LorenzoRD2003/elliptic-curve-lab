use num_bigint::BigUint;

use elliptic_algorithms_lab::elliptic_curves::ShortWeierstrassCurve;
use elliptic_algorithms_lab::fields::Fp17;
use elliptic_algorithms_lab::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};
use elliptic_algorithms_lab::visualization::{
    explain_graph_endomorphism_report, explain_isogeny_graph, explain_volcano_like_layers,
};

type F = Fp17;

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

    Ok(())
}
