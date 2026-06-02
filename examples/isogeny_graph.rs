use elliptic_algorithms_lab::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId, infer_volcano_like_layers,
};
use elliptic_algorithms_lab::{
    Field, Fp, ShortWeierstrassCurve, explain_isogeny_graph, explain_volcano_like_layers,
};

type F = Fp<17>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(1), F::from_i64(0))?;

    let graph = IsogenyGraphBuilder::new(curve, 2)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    println!("Sixth milestone: educational ℓ-isogeny graph explorer");
    println!("======================================================");
    println!();
    println!("{}", explain_isogeny_graph(&graph));
    println!();

    let report = graph.verify_locally()?;
    println!("Local verification report:");
    println!("{report:#?}");
    println!();

    let layers = infer_volcano_like_layers(&graph, IsogenyGraphNodeId(0));
    println!("{}", explain_volcano_like_layers(&graph, &layers));

    Ok(())
}
