use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

use crate::numerics::{
    cornacchia::{CornacchiaError, CornacchiaSolution},
    exact_square_root,
    hensel::{HenselLiftError, sqrt_mod_m},
};

use super::{root::*, validation::validate_coefficient_and_target};

/// Runs Cornacchia's algorithm for `x² + d y² = m` from a supplied root.
///
/// The supplied `r` must satisfy `r² ≡ -d mod m`. The algorithm normalizes it
/// to the representative `0 ≤ r ≤ m/2`, runs the partial Euclidean algorithm
/// until the current remainder `x` satisfies `x² < m`, and then tests whether
/// `(m - x²) / d` is a square.
///
/// Returns `Ok(Some(_))` when this root yields a solution, `Ok(None)` when the
/// final candidate is not an integer square, and `Err(_)` for invalid inputs or
/// an invalid supplied root.
///
/// Complexity: after validation, `Θ(log m)` Euclidean divisions plus one exact
/// integer-square-root test on an integer bounded by `m/d`.
pub fn cornacchia_with_root(
    d: &BigUint,
    m: &BigUint,
    r: &BigUint,
) -> Result<Option<CornacchiaSolution>, CornacchiaError> {
    validate_coefficient_and_target(d, m)?;
    validate_root(d, m, r)?;

    let mut previous = m.clone();
    let mut current = normalize_root(r, m);
    while &current * &current >= *m {
        let next = &previous % &current;
        previous = current;
        current = next;
    }

    let x_squared = &current * &current;
    let residual = m - &x_squared;
    if &residual % d != BigUint::zero() {
        return Ok(None);
    }

    let candidate = residual / d;
    let Some(y) = exact_square_root(&candidate) else {
        return Ok(None);
    };

    Ok(Some(CornacchiaSolution::new(current, y)))
}

/// Runs Cornacchia's algorithm against every square root of `-d` modulo `m`.
///
/// The returned solutions are distinct non-negative pairs `(x, y)` found by
/// running [`cornacchia_with_root`] on each modular square root of `-d`. For
/// non-square-free `m`, this is intentionally a candidate surface: it does not
/// promise to enumerate every integer solution of `x² + d y² = m`.
///
/// If `-d` has no square root modulo `m`, this returns an empty vector.
pub(crate) fn cornacchia_candidate_solutions(
    d: &BigUint,
    m: &BigUint,
) -> Result<Vec<CornacchiaSolution>, CornacchiaError> {
    validate_coefficient_and_target(d, m)?;

    let minus_d = -BigInt::from(d.clone());
    let roots = match sqrt_mod_m(&minus_d, m) {
        Ok(roots) => roots,
        Err(HenselLiftError::NoSquareRootModuloPrimePower) => return Ok(Vec::new()),
        Err(_) => return Err(CornacchiaError::ModularSquareRootFailure),
    };

    let mut solutions = Vec::new();
    for root in roots {
        if let Some(solution) = cornacchia_with_root(d, m, &root)? {
            solutions.push(solution);
        }
    }

    solutions.sort();
    solutions.dedup();
    Ok(solutions)
}

/// Returns the primitive Cornacchia candidate solutions for `x² + d y² = m`.
///
/// A candidate is primitive exactly when `gcd(x, y) = 1`.
pub fn cornacchia_primitive_solutions(
    d: &BigUint,
    m: &BigUint,
) -> Result<Vec<CornacchiaSolution>, CornacchiaError> {
    Ok(cornacchia_candidate_solutions(d, m)?
        .into_iter()
        .filter(CornacchiaSolution::is_primitive)
        .collect())
}
