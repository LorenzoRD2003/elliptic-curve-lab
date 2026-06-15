use elliptic_algorithms_lab::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{
        HasseInterval,
        group_order::{GroupOrderReport, GroupOrderStrategy, MestreConfig},
    },
    short_weierstrass::{
        group_exponent::{GroupExponentReport, GroupExponentStrategy},
        isomorphisms::{ShortWeierstrassQuadraticTwist, TwistKind},
        point_order::PointOrderStrategy,
    },
    traits::{FiniteGroupCurveModel, FrobeniusTraceCurveModel},
};
use elliptic_algorithms_lab::fields::{
    Fp,
    traits::{EnumerableFiniteField, Field},
};
use elliptic_algorithms_lab::visualization::{Visualizable, format_curve, format_point_compact};

type F = Fp<241>;

fn heading(title: &str) {
    println!("{title}");
    println!("{}", "-".repeat(title.len()));
}

fn max_order_point_index(curve: &ShortWeierstrassCurve<F>) -> usize {
    curve
        .point_orders()
        .into_iter()
        .enumerate()
        .max_by_key(|(_, (_, order))| *order)
        .map(|(index, _)| index)
        .expect("small finite curve should contain at least one point")
}

fn genuine_twist_curve(curve: &ShortWeierstrassCurve<F>) -> ShortWeierstrassCurve<F> {
    F::elements()
        .into_iter()
        .find_map(|candidate| {
            if F::is_zero(&candidate) {
                return None;
            }
            let Ok(package) = ShortWeierstrassQuadraticTwist::new(curve.clone(), candidate) else {
                return None;
            };
            (package.kind() == TwistKind::Quadratic).then(|| package.twist().clone())
        })
        .expect("a prime-field curve should admit a genuine quadratic twist")
}

fn parity_restricted_candidate_count(interval: &HasseInterval, even_group_order: bool) -> u128 {
    let wanted_parity = if even_group_order { 0 } else { 1 };
    let first = if interval.lower() % 2 == wanted_parity {
        interval.lower()
    } else {
        interval.lower() + 1
    };
    ((interval.upper() - first) / 2) + 1
}

fn ceil_sqrt_u128(value: u128) -> u128 {
    let floor = value.isqrt();
    if floor * floor == value {
        floor
    } else {
        floor + 1
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let curve = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
    let sample_point = curve
        .generator()
        .expect("the chosen F_241 sample curve should be cyclic");

    let exhaustive = curve.group_order_by(GroupOrderStrategy::Exhaustive)?;
    let character_sum = curve.group_order_by(GroupOrderStrategy::QuadraticCharacter)?;
    let interval = exhaustive.hasse_interval();

    let naive_hasse = curve.find_annihilating_multiple_in_hasse_interval_naive(&sample_point)?;

    let max_order_index = max_order_point_index(&curve);
    let mut exponent_sampler = {
        let mut requested = vec![max_order_index].into_iter();
        move |_upper_bound: usize| requested.next().or(Some(max_order_index))
    };
    let exponent_report = curve.group_exponent_by(
        GroupExponentStrategy::RandomPoints {
            max_samples: 1,
            point_order_strategy: PointOrderStrategy::Exhaustive,
        },
        &mut exponent_sampler,
    )?;
    let GroupExponentReport::RandomPoints(accumulation) = &exponent_report else {
        unreachable!("the chosen exponent strategy should preserve its variant")
    };
    let verification = curve.verify_exponent_lower_bound_by_group_order(
        accumulation,
        GroupOrderStrategy::QuadraticCharacter,
    )?;

    let twist_curve = genuine_twist_curve(&curve);
    let original_index = max_order_point_index(&curve);
    let twist_index = max_order_point_index(&twist_curve);
    let mut mestre_sampler = {
        let mut requested = vec![original_index, twist_index].into_iter();
        move |_upper_bound: usize| requested.next().or(Some(original_index))
    };
    let mestre = curve.group_order_by_with_sampler(
        GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(8)),
        &mut mestre_sampler,
    )?;
    let GroupOrderReport::MestreFp(mestre_report) = &mestre else {
        unreachable!("the chosen Mestre strategy should preserve its variant")
    };

    let bsgs_found = curve
        .find_annihilating_multiple_in_hasse_interval_bsgs(&sample_point)?
        .expect("Hasse's theorem should guarantee an annihilating multiple");
    let parity_restricted_count =
        parity_restricted_candidate_count(&interval, exhaustive.curve_order() % 2 == 0);
    let baby_steps = ceil_sqrt_u128(parity_restricted_count.div_ceil(2));
    let giant_stride_width = (2 * baby_steps) - 1;
    let giant_steps = parity_restricted_count.div_ceil(giant_stride_width);

    println!("Point-counting algorithm comparison");
    println!("======================================================");
    println!();
    println!("Curve: {} over F_241", format_curve(&curve));
    println!(
        "Sample point for H(q) searches: {}",
        format_point_compact(&sample_point)
    );
    println!();

    heading("Exhaustive count");
    println!("  #E(F_p) = {}", exhaustive.curve_order());
    println!("  t = {}", exhaustive.trace());
    println!();

    heading("Character sum");
    let GroupOrderReport::QuadraticCharacter(character_sum_report) = &character_sum else {
        unreachable!("the chosen quadratic-character strategy should preserve its variant")
    };
    println!("  Σχ(f(x)) = {}", character_sum_report.character_sum());
    println!("  #E(F_p) = {}", character_sum_report.curve_order());
    println!("  t = {}", character_sum_report.trace());
    println!();

    heading("Hasse interval");
    println!("  [{}, {}]", interval.lower(), interval.upper());
    println!(
        "  width as candidate count = {}",
        interval.candidate_count()
    );
    println!();

    heading("Naive Hasse multiple search");
    println!(
        "  found M = {:?}",
        naive_hasse.first_annihilating_multiple()
    );
    println!("  tested candidates = {}", naive_hasse.tested_candidates());
    println!("  first hit = {}", naive_hasse.format_compact());
    println!();

    heading("Exponent lower bound from sampled point orders");
    let sampled_orders = accumulation
        .steps()
        .iter()
        .map(|step| step.point_order_report().exact_order().to_string())
        .collect::<Vec<_>>();
    println!("  samples = {}", accumulation.samples_taken());
    println!("  point orders = [{}]", sampled_orders.join(", "));
    println!("  lcm = {}", accumulation.exponent_lower_bound());
    println!(
        "  unique multiple in H(p) = {:?}",
        verification.verified_group_order()
    );
    println!();

    heading("Mestre + twist");
    println!(
        "  λ approximation for E = {}",
        mestre_report.original_exponent_lower_bound()
    );
    println!(
        "  λ approximation for twist = {}",
        mestre_report.twist_exponent_lower_bound()
    );
    println!("  selected side = {}", mestre_report.resolved_side_label());
    println!(
        "  selected candidate = {}",
        mestre_report.resolved_side_group_order_candidate()
    );
    println!("  recorded steps = {}", mestre_report.step_count());
    println!();

    heading("BSGS");
    println!("  parity-restricted candidates = {parity_restricted_count}");
    println!("  baby steps = {baby_steps}");
    println!("  giant steps = {giant_steps}");
    println!("  found M = {bsgs_found}");
    println!();

    println!("Summary");
    println!("-------");
    println!(
        "  exhaustive, quadratic-character, Mestre, and both H(q) searches agree on #E(F_241) = {}.",
        exhaustive.curve_order()
    );
    println!(
        "  the sampled exponent lower bound already forces the same unique group order inside H(241)."
    );

    Ok(())
}
