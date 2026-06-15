use num_bigint::BigUint;
use std::collections::HashMap;
use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        hasse::search::config::{HasseBsgsConfig, HasseBsgsParity, HasseBsgsTraversal},
    },
    traits::{BigScalarGroupCurveModel, HasseIntervalSearchCurveModel},
};

pub(crate) fn find_annihilating_multiple_in_interval_bsgs_with_config_impl<
    C: HasseIntervalSearchCurveModel + ?Sized,
>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
    config: HasseBsgsConfig,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Eq + Hash,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    match config.known_parity() {
        HasseBsgsParity::Unknown => {
            let first = interval.lower();
            let candidate_offset = curve.mul_scalar_biguint(point, &BigUint::from(first))?;
            let target = curve.neg(&candidate_offset);
            let candidate_index = find_progression_solution_bsgs(
                curve,
                point,
                &target,
                interval.candidate_count(),
                config,
            )?;
            Ok(candidate_index.map(|index| {
                first
                    .checked_add(index)
                    .expect("Hasse candidate should stay in range")
            }))
        }
        known_parity => {
            let Some(first) = first_interval_candidate_with_parity(&interval, known_parity) else {
                return Ok(None);
            };
            let candidate_count = parity_restricted_candidate_count(&interval, first);
            let doubled_point = curve.double(point)?;
            let candidate_offset = curve.mul_scalar_biguint(point, &BigUint::from(first))?;
            let target = curve.neg(&candidate_offset);
            let candidate_index = find_progression_solution_bsgs(
                curve,
                &doubled_point,
                &target,
                candidate_count,
                config.with_known_parity(HasseBsgsParity::Unknown),
            )?;
            Ok(candidate_index.map(|index| {
                first
                    .checked_add(
                        index
                            .checked_mul(2)
                            .expect("parity index * 2 should stay in range"),
                    )
                    .expect("parity-restricted Hasse candidate should stay in range")
            }))
        }
    }
}

fn find_progression_solution_bsgs<C: HasseIntervalSearchCurveModel + ?Sized>(
    curve: &C,
    step_point: &C::Point,
    target: &C::Point,
    candidate_count: u128,
    config: HasseBsgsConfig,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Eq + Hash,
{
    if candidate_count == 0 {
        return Ok(None);
    }

    let r = choose_baby_step_count(candidate_count, config);
    let giant_stride_width = if config.uses_fast_negation() {
        r.checked_mul(2)
            .and_then(|double_r| double_r.checked_sub(1))
            .expect("2r - 1 should stay in range")
    } else {
        r
    };
    let mut baby_lookup = HashMap::with_capacity(r as usize);
    let mut baby = curve.identity();
    baby_lookup.insert(baby.clone(), 0u128);
    for j in 1..r {
        baby = curve.add(&baby, step_point)?;
        baby_lookup.entry(baby.clone()).or_insert(j);
    }

    let giant_stride = curve.mul_scalar_biguint(step_point, &BigUint::from(giant_stride_width))?;
    let context = BsgsSearchContext {
        step_point,
        target,
        candidate_count,
        r,
        giant_stride_width,
        giant_stride: &giant_stride,
        baby_lookup: &baby_lookup,
        config,
    };
    match config.traversal() {
        HasseBsgsTraversal::LeftToRight => {
            find_progression_solution_bsgs_left_to_right(curve, &context)
        }
        HasseBsgsTraversal::MiddleOut => find_progression_solution_bsgs_middle_out(curve, &context),
    }
}

struct BsgsSearchContext<'a, C: HasseIntervalSearchCurveModel + ?Sized>
where
    C::Point: Clone + Eq + Hash,
{
    step_point: &'a C::Point,
    target: &'a C::Point,
    candidate_count: u128,
    r: u128,
    giant_stride_width: u128,
    giant_stride: &'a C::Point,
    baby_lookup: &'a HashMap<C::Point, u128>,
    config: HasseBsgsConfig,
}

/// Baseline traversal: scan the giant-step blocks in increasing order.
///
/// This is the best fixed-instance baseline in the current codebase because
/// once the initial block is built, each further block costs exactly one
/// subtraction by the precomputed giant stride.
fn find_progression_solution_bsgs_left_to_right<C: HasseIntervalSearchCurveModel + ?Sized>(
    curve: &C,
    context: &BsgsSearchContext<'_, C>,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Eq + Hash,
{
    let block_count = context.candidate_count.div_ceil(context.giant_stride_width);
    let initial_candidate_index =
        initial_candidate_index_for_block(0, context.giant_stride_width, context.r, context.config);
    let initial_multiple =
        curve.mul_scalar_biguint(context.step_point, &BigUint::from(initial_candidate_index))?;
    let mut giant = curve.sub(context.target, &initial_multiple)?;

    for block_index in 0..block_count {
        if block_index > 0 {
            giant = curve.sub(&giant, context.giant_stride)?;
        }

        if let Some(candidate) = candidate_for_block_match(curve, context, &giant, block_index) {
            return Ok(Some(candidate));
        }
    }

    Ok(None)
}

/// Center-first traversal with two monotone frontiers.
///
/// The key implementation choice is that we do *not* alternate by
/// re-centering one giant-step state onto remote blocks. Instead we keep:
///
/// - one giant-step state expanding to the right from the center block
/// - one giant-step state expanding to the left from the center block
///
/// Therefore every newly visited block still costs only one group update, but
/// the visitation order is
///
/// `center, center + 1, center - 1, center + 2, center - 2, ...`
///
/// This is the operational form of the heuristic “look near `q + 1` first”.
fn find_progression_solution_bsgs_middle_out<C: HasseIntervalSearchCurveModel + ?Sized>(
    curve: &C,
    context: &BsgsSearchContext<'_, C>,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Eq + Hash,
{
    let block_count = context.candidate_count.div_ceil(context.giant_stride_width);
    let center_candidate_index = context.candidate_count.saturating_sub(1) / 2;
    let center_block_index = center_candidate_index / context.giant_stride_width;
    let initial_multiple = curve.mul_scalar_biguint(
        context.step_point,
        &BigUint::from(initial_candidate_index_for_block(
            center_block_index,
            context.giant_stride_width,
            context.r,
            context.config,
        )),
    )?;
    let center_giant = curve.sub(context.target, &initial_multiple)?;

    if let Some(candidate) =
        candidate_for_block_match(curve, context, &center_giant, center_block_index)
    {
        return Ok(Some(candidate));
    }

    let mut right_giant = center_giant.clone();
    let mut left_giant = center_giant;

    for offset in 1..block_count {
        if let Some(right_block_index) = center_block_index
            .checked_add(offset)
            .filter(|&index| index < block_count)
        {
            right_giant = curve.sub(&right_giant, context.giant_stride)?;
            if let Some(candidate) =
                candidate_for_block_match(curve, context, &right_giant, right_block_index)
            {
                return Ok(Some(candidate));
            }
        }

        if center_block_index >= offset {
            let left_block_index = center_block_index - offset;
            left_giant = curve.add(&left_giant, context.giant_stride)?;
            if let Some(candidate) =
                candidate_for_block_match(curve, context, &left_giant, left_block_index)
            {
                return Ok(Some(candidate));
            }
        }
    }

    Ok(None)
}

/// Interprets one giant-step state against the stored baby table.
///
/// For the fast-negation variant, the same giant-step state can certify two
/// candidate offsets because the lookup is tested against both `giant` and
/// `-giant`.
fn candidate_for_block_match<C: HasseIntervalSearchCurveModel + ?Sized>(
    curve: &C,
    context: &BsgsSearchContext<'_, C>,
    giant: &C::Point,
    block_index: u128,
) -> Option<u128>
where
    C::Point: Eq + Hash,
{
    if context.config.uses_fast_negation() {
        let giant_base = block_base(block_index, context.giant_stride_width);
        let giant_center = giant_base
            .checked_add(
                context
                    .r
                    .checked_sub(1)
                    .expect("fast-negation BSGS uses r >= 1"),
            )
            .expect("fast-negation giant-step center should stay in range");

        if let Some(&j) = context.baby_lookup.get(giant) {
            let candidate = giant_center
                .checked_add(j)
                .expect("center + j should stay in range");
            if candidate < context.candidate_count {
                return Some(candidate);
            }
        }

        if let Some(&j) = context.baby_lookup.get(&curve.neg(giant)) {
            let candidate = giant_center
                .checked_sub(j)
                .expect("center - j should stay in range");
            if candidate < context.candidate_count {
                return Some(candidate);
            }
        }
    } else if let Some(&j) = context.baby_lookup.get(giant) {
        let candidate = block_base(block_index, context.giant_stride_width)
            .checked_add(j)
            .expect("BSGS candidate should stay in range");
        if candidate < context.candidate_count {
            return Some(candidate);
        }
    }

    None
}

fn choose_baby_step_count(candidate_count: u128, config: HasseBsgsConfig) -> u128 {
    match config.traversal() {
        HasseBsgsTraversal::LeftToRight => {
            if config.uses_fast_negation() {
                ceil_sqrt_u128(candidate_count.div_ceil(2))
            } else {
                ceil_sqrt_u128(candidate_count)
            }
        }
        HasseBsgsTraversal::MiddleOut => {
            // This is an intentionally heuristic sizing rule: it aims the
            // baby-step count at the expected central distance rather than at
            // the full interval width. The benchmark suite keeps this choice
            // honest by comparing it both on one fixed instance and on a
            // center-heavy corpus.
            let expected_distance = expected_middle_distance(candidate_count);
            let target = if config.uses_fast_negation() {
                expected_distance / 2.0
            } else {
                expected_distance
            };
            ceil_sqrt_u128(target.ceil().max(1.0) as u128)
        }
    }
}

fn expected_middle_distance(candidate_count: u128) -> f64 {
    let alpha = 2.0 / (3.0 * std::f64::consts::PI);
    (candidate_count as f64) * alpha
}

fn block_base(block_index: u128, giant_stride_width: u128) -> u128 {
    block_index
        .checked_mul(giant_stride_width)
        .expect("block index * stride should stay in range")
}

fn initial_candidate_index_for_block(
    block_index: u128,
    giant_stride_width: u128,
    r: u128,
    config: HasseBsgsConfig,
) -> u128 {
    let base = block_base(block_index, giant_stride_width);
    if config.uses_fast_negation() {
        base.checked_add(r.checked_sub(1).expect("fast-negation BSGS uses r >= 1"))
            .expect("initial centered giant-step multiple should stay in range")
    } else {
        base
    }
}

fn first_interval_candidate_with_parity(
    interval: &HasseInterval,
    parity: HasseBsgsParity,
) -> Option<u128> {
    let wanted = match parity {
        HasseBsgsParity::Unknown => return Some(interval.lower()),
        HasseBsgsParity::Even => 0,
        HasseBsgsParity::Odd => 1,
    };
    let lower = interval.lower();
    let first = if lower % 2 == wanted {
        lower
    } else {
        lower.checked_add(1)?
    };
    interval.contains(first).then_some(first)
}

fn parity_restricted_candidate_count(interval: &HasseInterval, first: u128) -> u128 {
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
