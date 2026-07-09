use num_bigint::{BigInt, BigUint};
use num_traits::One;

use elliptic_algorithms_lab::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    endomorphisms::quadratic_orders::QuadraticDiscriminant,
    frobenius::cm::{CmTraceSignCurveModel, cm_absolute_trace_candidates},
    traits::{EnumerableCurveModel, FrobeniusTraceCurveModel, PointIndexSampler},
};
use elliptic_algorithms_lab::numerics::quadratic_forms::DiagonalBinaryQuadraticForm;
use elliptic_algorithms_lab::visualization::Visualizable;

type F29 = elliptic_algorithms_lab::fields::Fp29;

struct FixedIndexSampler {
    indices: Vec<usize>,
    cursor: usize,
}

impl FixedIndexSampler {
    fn new(indices: Vec<usize>) -> Self {
        Self { indices, cursor: 0 }
    }
}

impl PointIndexSampler for FixedIndexSampler {
    fn sample_index(&mut self, upper_bound: usize) -> Option<usize> {
        let index = *self.indices.get(self.cursor)?;
        self.cursor += 1;
        (index < upper_bound).then_some(index)
    }
}

fn find_distinguishing_point(
    curve: &ShortWeierstrassCurve<F29>,
    p: &BigUint,
    absolute_trace: &BigUint,
    sampler: &mut impl PointIndexSampler,
    max_attempts: usize,
) -> Result<Option<(BigInt, AffinePoint<F29>)>, Box<dyn std::error::Error>> {
    for _ in 0..max_attempts {
        let Some(point) = curve.random_point(sampler) else {
            return Ok(None);
        };

        if let Some(trace) =
            curve.cm_trace_from_absolute_trace_with_point(p, absolute_trace, &point)?
        {
            return Ok(Some((trace, point)));
        }
    }

    Ok(None)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Cornacchia and CM trace sign: Problem 2(e) miniature");
    println!("====================================================");
    println!();
    println!("We study E: y² = x³ - 35x - 98 and the CM relation 4p = a_p² + 7y².");
    println!("This example uses p = 29.");
    println!();

    let p = BigUint::from(29u8);
    let d = BigUint::from(7u8);
    let discriminant = QuadraticDiscriminant::new(BigInt::from(-7));
    let curve = ShortWeierstrassCurve::<F29>::new(F29::from_i64(-35), F29::from_i64(-98))?;

    println!("curve over F_29: {}", curve.format_compact());
    println!("-7 is a square modulo 29, so the CM trace should be non-zero.");
    println!();

    let form = DiagonalBinaryQuadraticForm::new(d)?;
    let representations = form.primitive_representations(&p)?;
    println!("Represent p by x² + 7y²:");
    for representation in &representations {
        println!("  29 = {}² + 7·{}²", representation.x(), representation.y());
    }
    println!();

    let trace_candidates = cm_absolute_trace_candidates(&discriminant, &p)?;
    println!("CM absolute trace candidates from 4p = a_p² + 7y²:");
    for candidate in &trace_candidates {
        println!(
            "  |a_p| = {}, v = {}",
            candidate.absolute_trace(),
            candidate.cm_multiplier()
        );
    }

    let Some(candidate) = trace_candidates.first() else {
        println!("No CM trace candidate was found.");
        return Ok(());
    };
    let absolute_trace = candidate.absolute_trace();
    println!();

    let mut sampler = FixedIndexSampler::new(vec![0, 3, 7, 11, 17, 23, 5, 1]);
    let Some((trace, witness)) =
        find_distinguishing_point(&curve, &p, absolute_trace, &mut sampler, 8)?
    else {
        println!("No sampled point distinguished the sign of a_p.");
        return Ok(());
    };

    let positive_order = &p + BigUint::one() - absolute_trace;
    let negative_order = &p + BigUint::one() + absolute_trace;
    let positive_kills = curve.cm_scalar_kills_point(&witness, &positive_order)?;
    let negative_kills = curve.cm_scalar_kills_point(&witness, &negative_order)?;

    println!("Determine the sign with a sampled point:");
    println!("  P = {}", witness.format_compact());
    println!("  [{positive_order}]P = O ? {positive_kills}");
    println!("  [{negative_order}]P = O ? {negative_kills}");
    println!("  therefore a_p = {trace}");
    println!();

    let enumerated_trace = curve.frobenius_trace()?;
    println!("Exhaustive check for this small example:");
    println!("  #E(F_29) = {}", enumerated_trace.curve_order());
    println!("  29 + 1 - #E(F_29) = {}", enumerated_trace.trace());

    Ok(())
}
