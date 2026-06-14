use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::HasseInterval;
use crate::elliptic_curves::order_from_multiple::mul_scalar_biguint;
use crate::elliptic_curves::traits::GroupCurveModel;
use num_bigint::BigUint;
use std::collections::HashMap;
use std::hash::Hash;

/// Traversal policy for the Hasse-interval baby-step/giant-step search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum HasseBsgsTraversal {
    LeftToRight,
    MiddleOut,
}

/// Optional parity information for the unknown annihilating multiple.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum HasseBsgsParity {
    Unknown,
    Even,
    Odd,
}

/// Internal configuration for Hasse-interval BSGS.
///
/// The current implementation uses only the default left-to-right traversal
/// and ignores the optimization toggles. The stored knobs are preparatory
/// scaffolding for future improvements from Lecture 7, §7.11.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct HasseBsgsConfig {
    pub traversal: HasseBsgsTraversal,
    pub use_fast_negation: bool,
    pub known_parity: HasseBsgsParity,
}

impl Default for HasseBsgsConfig {
    fn default() -> Self {
        Self {
            traversal: HasseBsgsTraversal::LeftToRight,
            use_fast_negation: true,
            known_parity: HasseBsgsParity::Unknown,
        }
    }
}

/// Searches one discrete Hasse interval with the baby-step/giant-step method
/// from Algorithm 7.9.
/// https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_notes7/
///
/// Given one point `P` and one explicit interval `H(q)`, this helper searches
/// for an integer `M ∈ H(q)` such that `[M]P = O`.
///
/// Complexity: Let `c = |H(q) ∩ ℤ|`. The current implementation chooses
/// `r = ceil(√c)` and `s = ceil(c/r)`, then performs:
///
/// - `Θ(r)` group additions to build the baby steps
/// - `Θ(1)` big-scalar multiplications to build `[a]P` and `[r]P`
/// - `Θ(s)` hash lookups and giant-step additions
///
/// Thus the dominant group-operation count is `Θ(r + s) = Θ(√c)`, which for
/// Hasse intervals is `Θ(∜q)`.
pub(crate) fn find_annihilating_multiple_in_interval_bsgs<C: GroupCurveModel + ?Sized>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Clone + Eq + Hash,
{
    find_annihilating_multiple_in_interval_bsgs_with_config(
        curve,
        point,
        interval,
        HasseBsgsConfig::default(),
    )
}

/// Internal configurable BSGS engine for one Hasse interval.
///
/// The current implementation preserves the existing left-to-right search
/// semantics regardless of configuration. Future work may use `config` to
/// enable middle-out traversal, fast-negation matching, or parity-aware
/// stepping without changing the call boundary.
pub(crate) fn find_annihilating_multiple_in_interval_bsgs_with_config<C: GroupCurveModel + ?Sized>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
    config: HasseBsgsConfig,
) -> Result<Option<u128>, CurveError>
where
    C::Point: Clone + Eq + Hash,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let _ = config;
    let candidate_count = interval.candidate_count();
    let r = if config.use_fast_negation {
        ceil_sqrt_u128(candidate_count.div_ceil(2))
    } else {
        ceil_sqrt_u128(candidate_count)
    };
    let giant_stride_width = if config.use_fast_negation {
        r.checked_mul(2)
            .and_then(|double_r| double_r.checked_sub(1))
            .expect("2r - 1 should stay in range")
    } else {
        r
    };
    let s = candidate_count.div_ceil(giant_stride_width);

    let mut baby_lookup = HashMap::with_capacity(r as usize);
    let mut baby = curve.identity();
    baby_lookup.insert(baby.clone(), 0u128);
    for j in 1..r {
        baby = curve.add(&baby, point)?;
        baby_lookup.entry(baby.clone()).or_insert(j);
    }

    let giant_stride = mul_scalar_biguint(curve, point, &BigUint::from(giant_stride_width))?;
    let initial_multiple = if config.use_fast_negation {
        interval
            .lower()
            .checked_add(r.checked_sub(1).expect("fast-negation BSGS uses r >= 1"))
            .expect("initial centered giant-step multiple should stay in range")
    } else {
        interval.lower()
    };
    let mut giant = mul_scalar_biguint(curve, point, &BigUint::from(initial_multiple))?;

    for i in 0..s {
        if config.use_fast_negation {
            let giant_base = interval
                .lower()
                .checked_add(
                    i.checked_mul(giant_stride_width)
                        .expect("i * (2r - 1) should stay in range"),
                )
                .expect("fast-negation giant-step base should stay in range");
            let giant_center = giant_base
                .checked_add(r.checked_sub(1).expect("fast-negation BSGS uses r >= 1"))
                .expect("fast-negation giant-step center should stay in range");

            if let Some(&j) = baby_lookup.get(&curve.neg(&giant)) {
                let candidate = giant_center
                    .checked_add(j)
                    .expect("center + j should stay in range");
                if interval.contains(candidate) {
                    return Ok(Some(candidate));
                }
            }

            if let Some(&j) = baby_lookup.get(&giant) {
                let candidate = giant_center
                    .checked_sub(j)
                    .expect("center - j should stay in range");
                if interval.contains(candidate) {
                    return Ok(Some(candidate));
                }
            }
        } else if let Some(&j) = baby_lookup.get(&curve.neg(&giant)) {
            let candidate = interval
                .lower()
                .checked_add(
                    i.checked_mul(giant_stride_width)
                        .expect("i * r should stay in range"),
                )
                .and_then(|base| base.checked_add(j))
                .expect("BSGS candidate inside the Hasse interval should stay in range");
            if interval.contains(candidate) {
                return Ok(Some(candidate));
            }
        }

        if i + 1 < s {
            giant = curve.add(&giant, &giant_stride)?;
        }
    }

    Ok(None)
}

fn ceil_sqrt_u128(value: u128) -> u128 {
    let floor = value.isqrt();
    if floor * floor == value {
        floor
    } else {
        floor + 1
    }
}
