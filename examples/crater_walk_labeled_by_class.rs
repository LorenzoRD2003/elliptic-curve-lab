use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};

use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::{
        binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    },
    traits::FrobeniusTraceCurveModel,
};
use elliptic_algorithms_lab::isogenies::{
    class_group_action::{ClassGroupActionPlan, CraterOrientationWitness},
    graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
};
use elliptic_algorithms_lab::visualization::Visualizable;

type F101 = elliptic_algorithms_lab::fields::Fp101;
type F37 = elliptic_algorithms_lab::fields::Fp37;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F101>::new(F101::from_i64(1), F101::from_i64(12))?;
    let frobenius_trace = curve.frobenius_trace()?;
    let frobenius_discriminant =
        frobenius_trace.trace() * frobenius_trace.trace() - BigInt::from(4 * 101);
    let graph = IsogenyGraphBuilder::new(curve.clone(), 3)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()?;

    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), BigUint::from(1u8))?;
    let ideal = PrimeNormIdeal::split(order, BigUint::from(3u8), BigUint::from(1u8))?;
    let ideal_norm = ideal.norm().clone();
    let ideal_root = ideal.root_mod_ell().clone();
    let order_discriminant = ideal.order().discriminant().value().clone();
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))?;
    let start = IsogenyGraphNodeId(0);

    let crater = graph.volcano_crater_report(ideal.norm())?;
    let labeled_walk = graph.labeled_crater_walk_report(&class_group, ideal, start)?;
    let orientation = CraterOrientationWitness::new(
        &crater,
        BTreeMap::from([
            (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
            (IsogenyGraphNodeId(1), IsogenyGraphNodeId(2)),
            (IsogenyGraphNodeId(2), IsogenyGraphNodeId(0)),
        ]),
    )?;
    let oriented_walk = labeled_walk.clone().with_user_orientation(orientation)?;
    let local_powers = [-1, 0, 1]
        .into_iter()
        .map(|exponent| {
            oriented_walk.apply_power_from(IsogenyGraphNodeId(0), BigInt::from(exponent))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let class_order_comparison =
        oriented_walk.compare_generator_order(&class_group, IsogenyGraphNodeId(0))?;
    let generated_subgroup =
        class_group.generated_subgroup(labeled_walk.form_label().reduced_form())?;
    let klein_curve = ShortWeierstrassCurve::<F37>::new(F37::from_i64(2), F37::from_i64(2))?;
    let klein_trace = klein_curve.frobenius_trace()?;
    let klein_frobenius_discriminant =
        klein_trace.trace() * klein_trace.trace() - BigInt::from(4 * 37);
    let klein_class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-84))?;
    let klein_order =
        ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-84), BigUint::from(1u8))?;
    let klein_first_ideal =
        PrimeNormIdeal::split(klein_order.clone(), BigUint::from(11u8), BigUint::from(2u8))?;
    let klein_second_ideal = PrimeNormIdeal::ramified(klein_order, BigUint::from(3u8))?;
    let klein_generators = [
        BinaryQuadraticForm::new(BigInt::from(2), BigInt::from(2), BigInt::from(11)),
        BinaryQuadraticForm::new(BigInt::from(3), BigInt::from(0), BigInt::from(7)),
    ];
    let klein_generated_subgroup =
        klein_class_group.generated_subgroup_by_set(&klein_generators)?;
    let klein_target = BinaryQuadraticForm::new(BigInt::from(5), BigInt::from(4), BigInt::from(5));
    let klein_action_plan = ClassGroupActionPlan::from_local_ideals(
        &klein_class_group,
        &[klein_first_ideal, klein_second_ideal],
        &klein_target,
    )?;

    println!("Crater walk labeled by an ideal/form class");
    println!("==========================================");
    println!();
    println!("Setup");
    println!("-----");
    println!("Curve: {} over F_101", curve.format_compact());
    println!("#E(F_101) = {}", frobenius_trace.curve_order());
    println!("Frobenius trace: t = {}", frobenius_trace.trace());
    println!("Frobenius discriminant: Δπ = {frobenius_discriminant}");
    println!("Local isogeny degree: ℓ = {ideal_norm}");
    println!("Quadratic order: discriminant D = {order_discriminant}");
    println!("Prime ideal: 𝔭 = ({ideal_norm}, ω - {ideal_root})");
    println!();

    println!("{}", crater.describe());
    println!();

    println!("{}", labeled_walk.describe());
    println!();
    println!("{}", oriented_walk.describe());
    println!();
    println!("Small local powers");
    println!("------------------");
    println!("These paths use the user-supplied crater orientation.");
    for power in &local_powers {
        println!("{}", power.format_compact());
    }
    println!();
    println!("{}", class_order_comparison.describe());
    println!();
    println!("{}", generated_subgroup.describe());
    println!();
    println!("What this certifies");
    println!("-------------------");
    println!("The ideal, the reduced form class, and the local crater prime are compatible.");
    println!("The recorded walk follows certified horizontal crater edges in graph order.");
    println!("The user-supplied orientation follows certified internal crater edges.");
    println!("The class-order comparison checks the observed oriented orbit length.");
    println!("The generated subgroup is algebraic; it is not yet a certified CM action.");
    println!();

    println!("Algebraic non-cyclic contrast");
    println!("-----------------------------");
    println!("Curve: {} over F_37", klein_curve.format_compact());
    println!("#E(F_37) = {}", klein_trace.curve_order());
    println!("Frobenius trace: t = {}", klein_trace.trace());
    println!("Frobenius discriminant: Δπ = {klein_frobenius_discriminant}");
    println!("Class-group discriminant: D = -84");
    println!();
    println!("{}", klein_generated_subgroup.describe());
    println!();
    println!("{}", klein_action_plan.describe());
    println!();
    println!("Interpretation");
    println!("--------------");
    println!("This Klein-style class group is generated by two independent local classes.");
    println!("A single oriented crater can only model one cyclic local direction.");
    println!(
        "A geometric multi-prime action would need one certified crater orientation per generator."
    );

    Ok(())
}
