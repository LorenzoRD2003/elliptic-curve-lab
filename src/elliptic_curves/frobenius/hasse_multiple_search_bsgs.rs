use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::HasseInterval;
use crate::elliptic_curves::order_from_multiple::mul_scalar_biguint;
use crate::elliptic_curves::traits::GroupCurveModel;
use num_bigint::BigUint;
use std::collections::HashMap;
use std::hash::Hash;

/// Searches one discrete Hasse interval with the baby-step/giant-step method
/// from Algorithm 7.9.
/// https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_notes7/
///
/// Given one point `P` and one explicit interval `H(q)`, this helper searches
/// for an integer `M ∈ H(q)` such that `[M]P = O`.
///
/// Complexity: Let `c = |H(q) ∩ Z|`. The current implementation chooses
/// `r = ceil(√c)` and `s = ceil(c/r)`, then performs:
///
/// - `Θ(r)` group additions to build the baby steps
/// - `Θ(1)` big-scalar multiplications to build `[a]P` and `[r]P`
/// - `Θ(s)` hash lookups and giant-step additions
///
/// Thus the dominant group-operation count is `Θ(r + s) = Θ(√c)`, which for
/// Hasse intervals is `Θ(∜q)`.
pub(crate) fn find_annihilating_multiple_in_interval_bsgs<C>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
) -> Result<Option<u128>, CurveError>
where
    C: GroupCurveModel + ?Sized,
    C::Point: Clone + Eq + Hash,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let candidate_count = interval.candidate_count();
    let r = ceil_sqrt_u128(candidate_count);
    let s = candidate_count.div_ceil(r);

    let mut baby_lookup = HashMap::with_capacity(r as usize);
    let mut baby = curve.identity();
    baby_lookup.insert(baby.clone(), 0u128);
    for j in 1..r {
        baby = curve.add(&baby, point)?;
        baby_lookup.entry(baby.clone()).or_insert(j);
    }

    let giant_stride = mul_scalar_biguint(curve, point, &BigUint::from(r))?;
    let mut giant = mul_scalar_biguint(curve, point, &BigUint::from(interval.lower()))?;

    for i in 0..s {
        if let Some(&j) = baby_lookup.get(&curve.neg(&giant)) {
            let candidate = interval
                .lower()
                .checked_add(i.checked_mul(r).expect("i * r should stay in range"))
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

#[cfg(test)]
mod tests {
    use super::find_annihilating_multiple_in_interval_bsgs;
    use crate::elliptic_curves::traits::HasseMultipleSearchCurveModel;
    use crate::elliptic_curves::{
        CurveModel, EnumerableCurveModel, GroupCurveModel, ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};

    type F241 = Fp<241>;

    #[test]
    fn bsgs_hasse_search_finds_an_annihilating_multiple_inside_the_same_hasse_interval() {
        let curve = ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
            .expect("valid curve");
        let point = curve
            .points()
            .into_iter()
            .find(|point| !curve.is_identity(point))
            .expect("small finite curve should contain a non-identity point");
        let interval =
            crate::elliptic_curves::HasseInterval::for_q(241).expect("valid Hasse interval");

        let naive = curve
            .find_annihilating_multiple_in_interval_naive(&point, interval.clone())
            .expect("naive Hasse search should succeed");
        let bsgs = find_annihilating_multiple_in_interval_bsgs(&curve, &point, interval)
            .expect("BSGS Hasse search should succeed")
            .expect("Hasse's theorem guarantees an annihilating multiple");

        assert!(naive.interval().contains(bsgs));
        assert!(curve.is_torsion_point(
            &point,
            u64::try_from(bsgs).expect("small-prime Hasse candidates fit in u64")
        ));
    }
}
