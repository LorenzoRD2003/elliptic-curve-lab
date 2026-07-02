use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

use crate::elliptic_curves::{
    CurveError,
    group_algorithms::cyclic_roots::{
        CyclicPrimeRootBezout, CyclicPrimeRootError, CyclicPrimeRootInput, CyclicPrimeRootOutcome,
        CyclicPrimeRootReport, CyclicPrimeRootStep, CyclicPrimeRootTrace,
    },
    traits::GroupCurveModel,
};
use crate::numerics::{extended_gcd_bigint, inverse_mod_biguint};

/// Computes one prime-degree root in a finite cyclic curve group.
///
/// The group is written additively. For a target `γ`, this attempts to solve
/// `[r]ρ = γ`, where `r` is prime and `|G| = a r^k`. If `k > 0`, the caller
/// supplies `δ`, a generator of the `r`-Sylow subgroup of order `r^k`, and the
/// algorithm is:
///
/// - compute `α = aγ` and `β = r^kγ`;
/// - find `x` by brute force with `α = xδ`;
/// - reject when `r` does not divide `x`;
/// - otherwise return `ρ = s(x/r)δ + tβ`, where `s a + t r^(k+1) = 1`.
///
/// If `k = 0`, multiplication by `r` is an automorphism of `G`, and the unique
/// root is `[r⁻¹ mod |G|]γ`.
pub(crate) fn compute_cyclic_prime_root_report<C: GroupCurveModel + ?Sized>(
    curve: &C,
    target: &C::Point,
    root_degree: BigUint,
    group_order: BigUint,
    sylow_generator: &C::Point,
) -> Result<CyclicPrimeRootReport<C::Point>, CyclicPrimeRootError>
where
    C::Point: Clone + PartialEq,
{
    if !curve.contains(target) || !curve.contains(sylow_generator) {
        return Err(CurveError::PointNotOnCurve.into());
    }

    let input = CyclicPrimeRootInput::from_group_order_and_prime(group_order, root_degree)?;
    validate_sylow_generator(curve, sylow_generator, &input)?;

    if !input.root_degree_divides_group_order() {
        let inverse = inverse_mod_biguint(input.root_degree(), input.prime_to_root_cofactor())
            .expect("prime r should be invertible modulo a when k = 0");
        let root = curve.mul_scalar(target, &inverse)?;
        let trace = CyclicPrimeRootTrace::new(None, None, None, None, Vec::new());
        return Ok(CyclicPrimeRootReport::new(
            input,
            target.clone(),
            sylow_generator.clone(),
            trace,
            CyclicPrimeRootOutcome::Root { root },
        ));
    }

    let alpha = curve.mul_scalar(target, input.prime_to_root_cofactor())?;
    let beta = curve.mul_scalar(target, input.sylow_order())?;
    let (discrete_log, steps) =
        brute_force_sylow_discrete_log(curve, sylow_generator, input.sylow_order(), &alpha)?;

    let outcome = if (&discrete_log % input.root_degree()).is_zero() {
        let bezout = bezout_for_root_formula(&input)?;
        let root = root_from_discrete_log_and_bezout(
            curve,
            sylow_generator,
            &beta,
            input.root_degree(),
            &discrete_log,
            &bezout,
        )?;
        let trace = CyclicPrimeRootTrace::new(
            Some(alpha),
            Some(beta),
            Some(discrete_log),
            Some(bezout),
            steps,
        );
        return Ok(CyclicPrimeRootReport::new(
            input,
            target.clone(),
            sylow_generator.clone(),
            trace,
            CyclicPrimeRootOutcome::Root { root },
        ));
    } else {
        CyclicPrimeRootOutcome::NoRoot
    };

    let trace = CyclicPrimeRootTrace::new(Some(alpha), Some(beta), Some(discrete_log), None, steps);
    Ok(CyclicPrimeRootReport::new(
        input,
        target.clone(),
        sylow_generator.clone(),
        trace,
        outcome,
    ))
}

fn validate_sylow_generator<C: GroupCurveModel + ?Sized>(
    curve: &C,
    sylow_generator: &C::Point,
    input: &CyclicPrimeRootInput,
) -> Result<(), CyclicPrimeRootError>
where
    C::Point: Clone + PartialEq,
{
    let killed_by_sylow_order = curve.mul_scalar(sylow_generator, input.sylow_order())?;
    if !curve.is_identity(&killed_by_sylow_order) {
        return Err(CyclicPrimeRootError::InvalidSylowGenerator {
            expected_order: input.sylow_order().clone(),
        });
    }

    if input.sylow_order() == &BigUint::one() {
        return Ok(());
    }

    let previous_sylow_order = input.sylow_order() / input.root_degree();
    let image = curve.mul_scalar(sylow_generator, previous_sylow_order)?;
    if curve.is_identity(&image) {
        return Err(CyclicPrimeRootError::InvalidSylowGenerator {
            expected_order: input.sylow_order().clone(),
        });
    }

    Ok(())
}

fn brute_force_sylow_discrete_log<C: GroupCurveModel + ?Sized>(
    curve: &C,
    sylow_generator: &C::Point,
    sylow_order: &BigUint,
    alpha: &C::Point,
) -> Result<(BigUint, Vec<CyclicPrimeRootStep<C::Point>>), CyclicPrimeRootError>
where
    C::Point: Clone + PartialEq,
{
    let mut candidate = BigUint::one();
    let mut candidate_multiple = sylow_generator.clone();
    let mut steps = Vec::new();

    while &candidate <= sylow_order {
        steps.push(CyclicPrimeRootStep::new(
            candidate.clone(),
            candidate_multiple.clone(),
        ));
        if &candidate_multiple == alpha {
            return Ok((candidate, steps));
        }

        candidate += BigUint::one();
        if &candidate <= sylow_order {
            candidate_multiple = curve.add(&candidate_multiple, sylow_generator)?;
        }
    }

    Err(CyclicPrimeRootError::MissingSylowDiscreteLog {
        sylow_order: sylow_order.clone(),
    })
}

fn bezout_for_root_formula(
    input: &CyclicPrimeRootInput,
) -> Result<CyclicPrimeRootBezout, CyclicPrimeRootError> {
    let next_sylow_order = input.sylow_order() * input.root_degree();
    let (gcd, s, t) = extended_gcd_bigint(
        BigInt::from(input.prime_to_root_cofactor().clone()),
        BigInt::from(next_sylow_order.clone()),
    );
    if gcd != BigInt::one() {
        return Err(CyclicPrimeRootError::MissingBezoutData {
            cofactor: input.prime_to_root_cofactor().clone(),
            next_sylow_order,
        });
    }
    Ok(CyclicPrimeRootBezout::new(
        s,
        t,
        input.prime_to_root_cofactor().clone(),
        next_sylow_order,
    ))
}

fn root_from_discrete_log_and_bezout<C: GroupCurveModel + ?Sized>(
    curve: &C,
    sylow_generator: &C::Point,
    beta: &C::Point,
    root_degree: &BigUint,
    discrete_log: &BigUint,
    bezout: &CyclicPrimeRootBezout,
) -> Result<C::Point, CyclicPrimeRootError>
where
    C::Point: Clone,
{
    let sylow_scalar = BigInt::from(discrete_log / root_degree) * bezout.s();
    let sylow_part = curve.mul_scalar_signed(sylow_generator, sylow_scalar)?;
    let beta_part = curve.mul_scalar_signed(beta, bezout.t().clone())?;
    Ok(curve.add(&sylow_part, &beta_part)?)
}
