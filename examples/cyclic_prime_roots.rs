//! Problem Set 2, Problem 3: prime-degree roots in a cyclic curve group.
//! https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_ps3/
//!
//! The exercise gives a cyclic group `G`, written additively, and asks for an
//! algorithm that solves `[r]ρ = γ` when `r` is prime. After writing
//! `|G| = a r^k`, the nontrivial case computes
//!
//! - `α = aγ`,
//! - `β = r^kγ`,
//! - a brute-force discrete logarithm `x` with `α = xδ`, where `δ` generates
//!   the `r`-Sylow subgroup,
//! - and then returns `ρ = s(x/r)δ + tβ`, with `s a + t r^(k+1) = 1`.
//!
//! This example instantiates the large curve from the statement over
//! `𝔽_p`, with `p = 2^255 - 19`. The original assignment chooses `c` from a
//! student ID; here we use the fixed value `c = 186 = 2·3·31` so that `cQ`
//! is visibly in the image of `[r]` for `r ∈ {2, 3, 31}`.

use std::str::FromStr;

use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    group_algorithms::cyclic_roots::{CyclicGroupPrimeRootCurveModel, CyclicPrimeRootOutcome},
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use elliptic_algorithms_lab::fields::{Fp25519, Fp25519Elem};
use num_bigint::BigUint;

type F = Fp25519;
type Point = AffinePoint<F>;
type Curve = ShortWeierstrassCurve<F>;

const COFACTOR_M: &str =
    "311269057089559665117126303786795451217418463436862985689835777395934466489";
const P_Y: &str = "3646051633135286488902046129458077014725501801396015176760137375427748642285";
const Q_X: &str = "43125933575059134974422288266359854378815207690220011740187158431378585841262";
const Q_Y: &str = "30438392960540783858586956956150489842875282144799753811252714114065692010946";

fn main() -> Result<(), CurveError> {
    let curve = problem_curve()?;
    let generator = curve.point(F::from_i64(99), fp(P_Y))?;
    let q = curve.point(fp(Q_X), fp(Q_Y))?;
    let group_order = problem_group_order();
    let c = BigUint::from(186u16);
    let target = curve.mul_scalar(&q, &c)?;

    assert!(
        curve.is_identity(&curve.mul_scalar(&generator, &group_order)?),
        "the supplied order n should annihilate the generator P"
    );

    println!("prime-degree roots on y^2 = x^3 + 31415926x + 27182818 over F_(2^255 - 19)");
    println!("fixed demonstration scalar: c = {c}");
    println!();
    println!("Target γ = cQ:");
    println_point("γ", &target);
    println!();

    for root_degree in [2u8, 3, 31] {
        run_prime_root_case(&curve, &generator, &target, &group_order, root_degree)?;
        println!();
    }

    Ok(())
}

fn problem_curve() -> Result<Curve, CurveError> {
    Curve::new(F::from_i64(31_415_926), F::from_i64(27_182_818))
}

fn problem_group_order() -> BigUint {
    BigUint::from(2u8) * BigUint::from(3u8) * BigUint::from(31u8) * bigint(COFACTOR_M)
}

fn run_prime_root_case(
    curve: &Curve,
    generator: &Point,
    target: &Point,
    group_order: &BigUint,
    root_degree: u8,
) -> Result<(), CurveError> {
    let r = BigUint::from(root_degree);
    let sylow_generator = curve.mul_scalar(generator, group_order / &r)?;
    let report = curve
        .cyclic_group_prime_root(target, r.clone(), group_order.clone(), &sylow_generator)
        .expect("the cyclic-root setup should be valid for the problem-set curve");

    println!("r = {root_degree}");
    println!("  |G| = a r^k with:");
    println!("    a   = {}", report.input().prime_to_root_cofactor());
    println!("    r^k = {}", report.input().sylow_order());
    println!("    k   = {}", report.input().sylow_exponent());
    println!("  δ = (n/r)P:");
    println_point("δ", &sylow_generator);

    match report.outcome() {
        CyclicPrimeRootOutcome::Root { root } => {
            let check = curve.mul_scalar(root, &r)?;
            assert_eq!(check, *target, "computed root must satisfy [r]R = γ");
            println!("  result: root found and verified by [r]R = γ");
            println_point("R", root);
        }
        CyclicPrimeRootOutcome::NoRoot => {
            println!("  result: no root");
            println!(
                "  certificate: x = {} is not divisible by r",
                report
                    .trace()
                    .discrete_log()
                    .expect("the nontrivial route records x before NoRoot")
            );
        }
    }

    Ok(())
}

fn fp(value: &str) -> Fp25519Elem {
    Fp25519Elem::from_biguint(&bigint(value))
}

fn bigint(value: &str) -> BigUint {
    BigUint::from_str(value).expect("decimal problem-set constant should parse")
}

fn println_point(name: &str, point: &Point) {
    match point {
        AffinePoint::Infinity => println!("  {name} = O"),
        AffinePoint::Finite { x, y } => {
            println!("  {name} = (");
            println!("    {},", x.to_biguint());
            println!("    {}", y.to_biguint());
            println!("  )");
        }
    }
}
