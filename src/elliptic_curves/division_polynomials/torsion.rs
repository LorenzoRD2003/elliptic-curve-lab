use crate::{
    elliptic_curves::{
        AffinePoint, CurveError, CurveModel, EnumerableCurveModel, GroupCurveModel,
        LiftXCoordinate, ShortWeierstrassCurve, point_has_exact_order, points_of_exact_order,
    },
    fields::{EnumerableFiniteField, SqrtField},
    polynomials::evaluation::evaluate_dense,
};

use super::{
    DivisionPolynomialError, DivisionPolynomialXCriterionKind,
    division_polynomial_x_criterion_kind, evaluate_division_polynomial_x_criterion,
    odd_division_polynomial,
};

/// Returns the rational `x`-coordinates in the base field that can correspond
/// to non-trivial `n`-torsion points when `n` is odd. This helper is intentionally
/// limited to small enumerable fields. It answers the question:
///
/// `Which rational x-coordinates can occur on points P with [n]P = O?`
///
/// Mathematical scope:
/// - it is defined only for odd `n`, because odd division polynomials lie in `F[x]`
/// - it returns only those roots `x in F` for which the curve actually has at
///   least one rational affine point above `x`
///
/// Algorithm:
/// 1. build the odd division polynomial `ψ_n(x)` once
/// 2. enumerate every `x in F`
/// 3. test whether `ψ_n(x) = 0`
/// 4. keep only those roots for which `x` lifts to at least one rational
///    affine point on the curve
///
/// Complexity:
/// - polynomial construction: same as [`odd_division_polynomial`], currently
///   about `Θ(n^5)` with naive dense multiplication
/// - root scan: one Horner evaluation per field element, so
///   `O(|F| · deg(ψ_n)) = O(|F| * n^2)`
///
/// If the crate later gains FFT-based multiplication,the overall cost
/// drops to roughly `Θ(n^3 log n + |F| * n^2)`.
fn rational_roots_of_odd_division_polynomial<F>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<F::Elem>, DivisionPolynomialError>
where
    F: EnumerableFiniteField + SqrtField,
{
    if n == 0 {
        return Err(DivisionPolynomialError::ZeroIndex);
    } else if n.is_multiple_of(2) {
        return Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n });
    }

    let polynomial = odd_division_polynomial(curve, n)?;
    let mut roots = Vec::new();
    for x in F::elements() {
        let value = evaluate_dense(&polynomial, &x)?;
        if F::is_zero(&value) && curve.point_from_x(x.clone()).is_some() {
            roots.push(x);
        }
    }
    Ok(roots)
}

/// Returns the rational affine points found by lifting the rational roots of
/// the odd division polynomial `ψ_n(x)`. It answers the question
///
/// `Which rational points can we recover by solving ψ_n(x) = 0 and then
/// lifting those x-coordinates back to the curve equation?`
///
/// Current implementation strategy:
///
/// 1. compute the odd division polynomial `ψ_n(x)`
/// 2. enumerate every `x in F` and keep the rational roots
/// 3. for each such `x`, enumerate every `y in F`
/// 4. retain the affine points `(x, y)` that satisfy
///    `y^2 = x^3 + ax + b`
///
/// This implementation chooses the fully enumerable path, so it currently
/// requires `F` to expose its full element set honestly. That makes the helper
/// suitable for the same small finite educational settings as the rest of the
/// milestone-7 torsion tooling.
///
/// Mathematical scope:
///
/// - this helper is intentionally restricted to odd `n`, because only then
///   does `ψ_n` live directly in `F[x]`
/// - even indices are rejected honestly via
///   [`DivisionPolynomialError::EvenIndexRequiresYFactor`]
///
/// Complexity:
///
/// - constructing `ψ_n`: about `Θ(n^5)` with the current naive dense backend
/// - scanning `x`: `O(|F| · n^2)` via Horner evaluation
/// - lifting each surviving `x` by brute-force `y` search: `O(r · |F|)`,
///   where `r` is the number of rational `x`-roots
///
/// So the full cost is roughly
/// `Θ(n^5 + |F| · n^2 + r · |F|)`.
fn torsion_candidates_from_odd_division_polynomial<F>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError>
where
    F: EnumerableFiniteField,
{
    if n == 0 {
        return Err(DivisionPolynomialError::ZeroIndex);
    }

    if n.is_multiple_of(2) {
        return Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n });
    }

    let polynomial = odd_division_polynomial(curve, n)?;
    let mut points = Vec::new();

    for x in F::elements() {
        let value = evaluate_dense(&polynomial, &x)?;
        if !F::is_zero(&value) {
            continue;
        }

        for y in F::elements() {
            let candidate = AffinePoint::<F>::new(x.clone(), y);
            if curve.contains(&candidate) {
                points.push(candidate);
            }
        }
    }

    Ok(points)
}

/// Returns the rational `x`-coordinates in the base field that can correspond
/// to affine points `P` with `ψ_n(P) = 0` when `n` is even.
///
/// For even indices the division polynomial has the form `ψ_n = y ε_n(x)`
/// with `ε_n(x) ∈ F[x]`. Therefore the vanishing condition at a rational point
/// `P = (x, y)` is `y = 0`, or `ε_n(x) = 0`.
///
/// This is why the function returns *candidates* rather than claiming exact
/// `n`-torsion from the polynomial alone.
///
/// Algorithm:
/// 1. enumerate every rational affine point of the curve by scanning
///    `x, y ∈ F`
/// 2. for each point, test whether `y = 0` or `ε_n(x) = 0`
/// 3. keep the distinct `x`-coordinates of the points that satisfy that
///    condition
fn rational_x_candidates_from_even_division_polynomial<F: EnumerableFiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<F::Elem>, DivisionPolynomialError> {
    if n == 0 {
        return Err(DivisionPolynomialError::ZeroIndex);
    }
    if !n.is_multiple_of(2) {
        return Err(DivisionPolynomialError::UnsupportedIndex { n });
    }
    let mut xs = Vec::new();

    for x in F::elements() {
        let even_factor_zero = evaluate_division_polynomial_x_criterion(curve, n, &x)
            .is_ok_and(|value| F::is_zero(&value));

        for y in F::elements() {
            let point = AffinePoint::<F>::new(x.clone(), y);
            if !curve.contains(&point) {
                continue;
            }
            let is_y_zero = match &point {
                AffinePoint::Finite { y, .. } => F::is_zero(y),
                AffinePoint::Infinity => false,
            };
            if (is_y_zero || even_factor_zero) && !xs.iter().any(|existing| F::eq(existing, &x)) {
                xs.push(x.clone());
            }
        }
    }
    Ok(xs)
}

/// Returns the rational affine points `P` with `ψ_n(P) = 0` when `n` is even.
///
/// For even `n`, the curve-side vanishing condition is
/// `ψ_n(P) = y(P) ε_n(x(P)) = 0`.
///
/// So a rational affine point is returned exactly when either:
/// - its `y`-coordinate is `0`, or
/// - its `x`-coordinate is a root of the even division-polynomial factor
///   `ε_n(x)`
fn torsion_candidates_from_even_division_polynomial<F: EnumerableFiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
    if n == 0 {
        return Err(DivisionPolynomialError::ZeroIndex);
    }
    if !n.is_multiple_of(2) {
        return Err(DivisionPolynomialError::UnsupportedIndex { n });
    }
    let mut points = Vec::new();

    for x in F::elements() {
        let even_factor_zero = evaluate_division_polynomial_x_criterion(curve, n, &x)
            .is_ok_and(|value| F::is_zero(&value));

        for y in F::elements() {
            let point = AffinePoint::<F>::new(x.clone(), y);
            if !curve.contains(&point) {
                continue;
            }
            let is_y_zero = match &point {
                AffinePoint::Finite { y, .. } => F::is_zero(y),
                AffinePoint::Infinity => false,
            };
            if is_y_zero || even_factor_zero {
                points.push(point);
            }
        }
    }
    Ok(points)
}

/// Returns the rational `x`-coordinates in the base field that can correspond
/// to rational affine points annihilated by the division polynomial `ψ_n`.
///
/// This is the main public `x`-coordinate query for the milestone-7 torsion
/// layer. It dispatches by parity:
///
/// - if `n` is odd, it solves `ψ_n(x) = 0` in `F`
/// - if `n` is even, it solves the curve-side condition `y ε_n(x) = 0`, so `x` is
///   accepted when it appears on a rational point with `y = 0` or with `ε_n(x) = 0`
pub fn rational_x_candidates_for_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<F::Elem>, DivisionPolynomialError> {
    match division_polynomial_x_criterion_kind(n)? {
        DivisionPolynomialXCriterionKind::OddDivisionPolynomial => {
            rational_roots_of_odd_division_polynomial(curve, n)
        }
        DivisionPolynomialXCriterionKind::EvenYStrippedFactor => {
            rational_x_candidates_from_even_division_polynomial(curve, n)
        }
    }
}

/// Returns the rational affine points found from the division-polynomial
/// vanishing condition `ψ_n(P) = 0`.
///
/// - for odd `n`, it lifts the rational roots of `ψ_n(x)`
/// - for even `n`, it enumerates the rational affine points satisfying
///   `y(P) = 0` or `ε_n(x(P)) = 0`
///
/// These are *candidates* rather than a certified exact-order set.
/// In particular, for even `n` the `y = 0` branch can contribute points coming
/// from lower-order `2`-torsion.
pub fn torsion_candidates_from_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
    match division_polynomial_x_criterion_kind(n)? {
        DivisionPolynomialXCriterionKind::OddDivisionPolynomial => {
            torsion_candidates_from_odd_division_polynomial(curve, n)
        }
        DivisionPolynomialXCriterionKind::EvenYStrippedFactor => {
            torsion_candidates_from_even_division_polynomial(curve, n)
        }
    }
}

/// Returns the rational affine points kept after a light torsion validation
/// pass on top of the division-polynomial candidates.
///
/// It differs from [`torsion_candidates_from_division_polynomial`] only in the
/// even case:
/// - odd `n`: every candidate is returned unchanged
/// - even `n`: candidates with `y = 0` are retained only if they are actually
///   killed by `[n]`
///
/// Educational note: on a short-Weierstrass curve over characteristic different from
/// `2`, any affine point with `y = 0` is automatically `2`-torsion, so for even `n`
/// such a point is indeed `n`-torsion. The explicit validation step is kept
/// anyway because it makes the branch logic visible and documents the intended
/// semantic difference between “candidate from `ψ_n(P)=0`” and “point we are
/// willing to present as torsion-confirmed”.
pub fn torsion_points_from_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
    let candidates = torsion_candidates_from_division_polynomial(curve, n)?;
    if !n.is_multiple_of(2) {
        return Ok(candidates);
    }
    let mut points = Vec::new();
    for point in candidates {
        match &point {
            AffinePoint::Infinity => {}
            AffinePoint::Finite { y, .. } => {
                if !F::is_zero(y) || curve.is_torsion_point(&point, n as u64) {
                    points.push(point);
                }
            }
        }
    }
    Ok(points)
}

/// Returns the rational affine points of exact order `n` recovered from the
/// division-polynomial search pipeline.
///
/// This helper first computes [`torsion_points_from_division_polynomial`] and
/// then validates exact order `n` point-by-point with the generic torsion
/// checker from `elliptic_curves::torsion`.
pub fn exact_n_torsion_points_from_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
    let candidates = torsion_points_from_division_polynomial(curve, n)?;
    let mut points = Vec::new();

    for point in candidates {
        if point_has_exact_order(curve, &point, n).map_err(|error| match error {
            CurveError::InvalidTorsionOrder { order: _ } => DivisionPolynomialError::ZeroIndex,
            other => DivisionPolynomialError::Curve(other),
        })? {
            points.push(point);
        }
    }

    Ok(points)
}

/// Comparison summary between division-polynomial-based torsion recovery and
/// direct small-group enumeration.
///
/// Interpretation:
///
/// - `polynomial_candidate_count` counts the affine points returned by the raw
///   division-polynomial candidate pipeline
/// - `polynomial_n_torsion_count` counts the points kept after the public
///   torsion-validation layer
/// - `enumerated_n_torsion_count` counts all non-identity affine points
///   satisfying `[n]P = O` by direct group traversal
/// - `exact_order_polynomial_count` counts the exact-order points recovered
///   from the division-polynomial pipeline
/// - `exact_order_enumerated_count` counts the exact-order points found by
///   exhaustive group enumeration
/// - `missing_from_polynomial` and `extra_from_polynomial` compare the
///   exact-order sets
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TorsionComparisonReport<F: EnumerableFiniteField + SqrtField> {
    pub n: usize,
    pub polynomial_candidate_count: usize,
    pub polynomial_n_torsion_count: usize,
    pub enumerated_n_torsion_count: usize,
    pub exact_order_polynomial_count: usize,
    pub exact_order_enumerated_count: usize,
    pub missing_from_polynomial: Vec<AffinePoint<F>>,
    pub extra_from_polynomial: Vec<AffinePoint<F>>,
}

/// Compares division-polynomial torsion recovery against direct exhaustive
/// group enumeration.
///
/// This helper is meant for the current small finite educational setting where
/// the full rational point set can be enumerated honestly. It summarizes:
///
/// - how many affine points arise as raw division-polynomial candidates
/// - how many survive the public torsion validation layer
/// - how many are found by direct `[n]P = O` enumeration
/// - how the exact-order-`n` sets compare
///
/// The `missing_from_polynomial` and `extra_from_polynomial` fields compare
/// the exact-order sets, since that is the most informative notion of
/// correctness once candidates and lower-order contamination are taken into
/// account.
pub fn compare_division_polynomial_torsion_with_enumeration<
    F: EnumerableFiniteField + SqrtField,
>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<TorsionComparisonReport<F>, DivisionPolynomialError> {
    let polynomial_candidates = torsion_candidates_from_division_polynomial(curve, n)?;
    let polynomial_n_torsion = torsion_points_from_division_polynomial(curve, n)?;
    let exact_order_polynomial = exact_n_torsion_points_from_division_polynomial(curve, n)?;

    let enumerated_n_torsion: Vec<_> = curve
        .points()
        .into_iter()
        .filter(|point| !curve.is_identity(point) && curve.is_torsion_point(point, n as u64))
        .collect();

    let exact_order_enumerated =
        points_of_exact_order(curve, n).map_err(DivisionPolynomialError::Curve)?;

    let mut missing_from_polynomial = Vec::new();
    for point in &exact_order_enumerated {
        if !exact_order_polynomial
            .iter()
            .any(|candidate| candidate == point)
        {
            missing_from_polynomial.push(point.clone());
        }
    }

    let mut extra_from_polynomial = Vec::new();
    for point in &exact_order_polynomial {
        if !exact_order_enumerated
            .iter()
            .any(|candidate| candidate == point)
        {
            extra_from_polynomial.push(point.clone());
        }
    }

    Ok(TorsionComparisonReport {
        n,
        polynomial_candidate_count: polynomial_candidates.len(),
        polynomial_n_torsion_count: polynomial_n_torsion.len(),
        enumerated_n_torsion_count: enumerated_n_torsion.len(),
        exact_order_polynomial_count: exact_order_polynomial.len(),
        exact_order_enumerated_count: exact_order_enumerated.len(),
        missing_from_polynomial,
        extra_from_polynomial,
    })
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        EnumerableCurveModel, FiniteGroupCurveModel,
        elliptic_curves::division_polynomials::evaluate_division_polynomial_at_point,
        fields::{EnumerableFiniteField, Field, Fp, SqrtField},
        proptest_support::non_singular_short_weierstrass_curve,
    };

    type F17 = Fp<17>;
    type F23 = Fp<23>;

    fn unique_x_coordinates_of_rational_n_torsion<F: EnumerableFiniteField + SqrtField>(
        curve: &ShortWeierstrassCurve<F>,
        n: usize,
    ) -> Vec<F::Elem> {
        let mut xs = Vec::new();
        for point in curve.points() {
            let AffinePoint::Finite { ref x, .. } = point else {
                continue;
            };
            if !curve.is_torsion_point(&point, n as u64) {
                continue;
            }
            if xs.iter().any(|existing| F::eq(existing, x)) {
                continue;
            }
            xs.push(x.clone());
        }
        xs
    }

    fn same_x_set<F: Field>(left: &[F::Elem], right: &[F::Elem]) -> bool {
        left.len() == right.len()
            && left.iter().all(|x| right.iter().any(|y| F::eq(x, y)))
            && right.iter().all(|x| left.iter().any(|y| F::eq(x, y)))
    }

    #[test]
    fn rational_roots_of_psi_3_match_x_coordinates_of_three_torsion() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let roots =
            rational_roots_of_odd_division_polynomial(&curve, 3).expect("psi_3 roots should exist");
        let expected = unique_x_coordinates_of_rational_n_torsion(&curve, 3);

        assert!(same_x_set::<F23>(&roots, &expected));
    }

    #[test]
    fn roots_lift_to_curve_points_when_rhs_is_square() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let roots =
            rational_roots_of_odd_division_polynomial(&curve, 3).expect("psi_3 roots should exist");
        let lifted =
            torsion_candidates_from_division_polynomial(&curve, 3).expect("psi_3 lift should work");
        let x = F23::from_i64(8);

        assert!(roots.iter().any(|root| F23::eq(root, &x)));
        assert!(curve.point_from_x(x).is_some());
        assert!(lifted.iter().all(|point| curve.contains(point)));
        assert!(
            lifted
                .iter()
                .all(|point| matches!(point, AffinePoint::Finite { .. }))
        );
        assert!(lifted.iter().any(
            |point| matches!(point, AffinePoint::Finite { x: root, .. } if F23::eq(root, &x))
        ));
    }

    #[test]
    fn non_liftable_roots_are_not_returned_as_rational_points() {
        type F5 = Fp<5>;

        let curve = ShortWeierstrassCurve::<F5>::new(F5::zero(), F5::one())
            .expect("curve should be non-singular");
        let polynomial =
            odd_division_polynomial(&curve, 3).expect("psi_3 should be available as a polynomial");
        let non_liftable_root = F5::one();
        let rhs = curve.rhs(&non_liftable_root);
        let rational_roots =
            rational_roots_of_odd_division_polynomial(&curve, 3).expect("psi_3 roots should exist");
        let lifted =
            torsion_candidates_from_division_polynomial(&curve, 3).expect("psi_3 lift should work");

        assert!(F5::is_zero(
            &evaluate_dense(&polynomial, &non_liftable_root)
                .expect("raw psi_3 root should evaluate")
        ));
        assert!(curve.point_from_x(non_liftable_root).is_none());
        assert!(curve.points_from_x(non_liftable_root).is_none());
        assert_eq!(rhs, F5::from_i64(2));
        assert!(
            !rational_roots
                .iter()
                .any(|root| F5::eq(root, &non_liftable_root))
        );
        assert!(lifted.iter().all(|point| match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, .. } => !F5::eq(x, &non_liftable_root),
        }));
    }

    #[test]
    fn public_x_candidates_dispatch_by_parity() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let odd = rational_x_candidates_for_division_polynomial(&curve, 3)
            .expect("odd x-candidates should exist");
        let odd_expected =
            rational_roots_of_odd_division_polynomial(&curve, 3).expect("psi_3 roots should exist");
        let even = rational_x_candidates_for_division_polynomial(&curve, 4)
            .expect("even x-candidates should exist");
        let even_expected = rational_x_candidates_from_even_division_polynomial(&curve, 4)
            .expect("psi_4 x-candidates should exist");

        assert!(same_x_set::<F23>(&odd, &odd_expected));
        assert!(same_x_set::<F23>(&even, &even_expected));
    }

    #[test]
    fn rational_roots_match_recursive_nine_torsion_x_coordinates_over_f17() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(1), F17::from_i64(1))
            .expect("curve should be non-singular");

        let roots =
            rational_roots_of_odd_division_polynomial(&curve, 9).expect("psi_9 roots should exist");
        let expected = unique_x_coordinates_of_rational_n_torsion(&curve, 9);

        assert!(same_x_set::<F17>(&roots, &expected));
    }

    #[test]
    fn rational_roots_can_be_empty_when_no_rational_odd_torsion_exists() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let roots =
            rational_roots_of_odd_division_polynomial(&curve, 5).expect("psi_5 roots should exist");

        assert!(roots.is_empty());
    }

    #[test]
    fn rational_roots_reject_even_and_zero_indices() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            rational_roots_of_odd_division_polynomial(&curve, 0),
            Err(DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            rational_roots_of_odd_division_polynomial(&curve, 6),
            Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n: 6 })
        );
    }

    fn same_point_set<F: Field>(left: &[AffinePoint<F>], right: &[AffinePoint<F>]) -> bool {
        left.len() == right.len()
            && left
                .iter()
                .all(|point| right.iter().any(|other| point == other))
            && right
                .iter()
                .all(|point| left.iter().any(|other| point == other))
    }

    #[test]
    fn polynomial_three_torsion_matches_enumerated_three_torsion() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let candidates =
            torsion_candidates_from_division_polynomial(&curve, 3).expect("psi_3 lift should work");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter(|point| !curve.is_identity(point) && curve.is_torsion_point(point, 3))
            .collect();

        assert!(same_point_set::<F23>(&candidates, &expected));
    }

    #[test]
    fn torsion_candidates_match_recursive_nine_torsion_points_over_f17() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(1), F17::from_i64(1))
            .expect("curve should be non-singular");

        let candidates =
            torsion_candidates_from_division_polynomial(&curve, 9).expect("psi_9 lift should work");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter(|point| !curve.is_identity(point) && curve.is_torsion_point(point, 9))
            .collect();

        assert!(same_point_set::<F17>(&candidates, &expected));
    }

    #[test]
    fn public_torsion_candidates_dispatch_by_parity() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let odd = torsion_candidates_from_division_polynomial(&curve, 3)
            .expect("odd torsion candidates should exist");
        let odd_expected = torsion_candidates_from_odd_division_polynomial(&curve, 3)
            .expect("odd torsion candidates should exist");
        let even = torsion_candidates_from_division_polynomial(&curve, 6)
            .expect("even torsion candidates should exist");
        let even_expected = torsion_candidates_from_even_division_polynomial(&curve, 6)
            .expect("even torsion candidates should exist");

        assert!(same_point_set::<F23>(&odd, &odd_expected));
        assert!(same_point_set::<F23>(&even, &even_expected));
    }

    #[test]
    fn torsion_candidates_public_api_rejects_only_zero_index() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            torsion_candidates_from_division_polynomial(&curve, 0),
            Err(DivisionPolynomialError::ZeroIndex)
        );
    }

    #[test]
    fn even_x_candidates_match_points_where_psi_n_vanishes_over_f23() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let xs = rational_x_candidates_from_even_division_polynomial(&curve, 4)
            .expect("psi_4 x-candidates should exist");
        let mut expected = Vec::new();

        for point in curve.points() {
            let AffinePoint::Finite { x, .. } = point else {
                continue;
            };

            let value = evaluate_division_polynomial_at_point(&curve, 4, &point)
                .expect("psi_4(P) should evaluate");
            if F23::is_zero(&value) && !expected.iter().any(|existing| F23::eq(existing, &x)) {
                expected.push(x);
            }
        }

        assert!(same_x_set::<F23>(&xs, &expected));
    }

    #[test]
    fn even_torsion_candidates_match_points_where_psi_n_vanishes_over_f23() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let candidates = torsion_candidates_from_even_division_polynomial(&curve, 6)
            .expect("psi_6 candidates should exist");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter_map(|point| match point {
                AffinePoint::Infinity => None,
                finite => {
                    let value = evaluate_division_polynomial_at_point(&curve, 6, &finite)
                        .expect("psi_6(P) should evaluate");
                    F23::is_zero(&value).then_some(finite)
                }
            })
            .collect();

        assert!(same_point_set::<F23>(&candidates, &expected));
    }

    #[test]
    fn torsion_points_validate_even_y_zero_branch() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let candidates = torsion_candidates_from_division_polynomial(&curve, 6)
            .expect("psi_6 candidates should exist");
        let points = torsion_points_from_division_polynomial(&curve, 6)
            .expect("psi_6 torsion points should exist");

        let expected: Vec<_> = candidates
            .iter()
            .filter_map(|point| match point {
                AffinePoint::Infinity => None,
                AffinePoint::Finite { y, .. } => {
                    if !F23::is_zero(y) || curve.is_torsion_point(point, 6) {
                        Some(point.clone())
                    } else {
                        None
                    }
                }
            })
            .collect();

        assert!(same_point_set::<F23>(&points, &expected));
    }

    #[test]
    fn exact_order_three_points_exclude_identity() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");
        let identity = curve.identity();

        let exact_points = exact_n_torsion_points_from_division_polynomial(&curve, 3)
            .expect("exact order-three torsion should exist");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter(|point| !curve.is_identity(point) && curve.point_order(point) == Some(3))
            .collect();

        assert!(!exact_points.contains(&identity));
        assert!(same_point_set::<F23>(&exact_points, &expected));
    }

    #[test]
    fn n_torsion_and_exact_order_n_are_distinguished() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let candidates = torsion_points_from_division_polynomial(&curve, 12)
            .expect("torsion candidates should exist");
        let exact_points = exact_n_torsion_points_from_division_polynomial(&curve, 12)
            .expect("exact order-twelve torsion should exist");
        let expected: Vec<_> = curve
            .points()
            .into_iter()
            .filter(|point| !curve.is_identity(point) && curve.point_order(point) == Some(12))
            .collect();

        assert!(candidates.len() >= exact_points.len());
        assert!(same_point_set::<F23>(&exact_points, &expected));
        assert!(
            candidates
                .iter()
                .any(|point| curve.point_order(point) == Some(2))
        );
    }

    #[test]
    fn m6_torsion_enumeration_and_m7_polynomial_torsion_agree_on_small_curve() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        let report = compare_division_polynomial_torsion_with_enumeration(&curve, 3)
            .expect("comparison report should build");

        assert_eq!(report.polynomial_candidate_count, 2);
        assert_eq!(
            report.polynomial_n_torsion_count,
            report.enumerated_n_torsion_count
        );
        assert_eq!(
            report.exact_order_polynomial_count,
            report.exact_order_enumerated_count
        );
        assert!(report.missing_from_polynomial.is_empty());
        assert!(report.extra_from_polynomial.is_empty());
    }

    #[test]
    fn exact_n_torsion_points_reject_zero_index() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            exact_n_torsion_points_from_division_polynomial(&curve, 0),
            Err(DivisionPolynomialError::ZeroIndex)
        );
    }

    #[test]
    fn torsion_points_match_candidates_for_odd_indices() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(1), F17::from_i64(1))
            .expect("curve should be non-singular");

        let candidates = torsion_candidates_from_division_polynomial(&curve, 9)
            .expect("psi_9 candidates should exist");
        let points = torsion_points_from_division_polynomial(&curve, 9)
            .expect("psi_9 torsion points should exist");

        assert!(same_point_set::<F17>(&points, &candidates));
    }

    #[test]
    fn even_candidate_helpers_reject_zero_and_odd_indices() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            rational_x_candidates_from_even_division_polynomial(&curve, 0),
            Err(DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            rational_x_candidates_from_even_division_polynomial(&curve, 3),
            Err(DivisionPolynomialError::UnsupportedIndex { n: 3 })
        );
        assert_eq!(
            torsion_candidates_from_even_division_polynomial(&curve, 0),
            Err(DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            torsion_candidates_from_even_division_polynomial(&curve, 3),
            Err(DivisionPolynomialError::UnsupportedIndex { n: 3 })
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn property_exact_odd_torsion_matches_enumeration(
            curve in non_singular_short_weierstrass_curve::<17>(),
        ) {
            let polynomial_points = exact_n_torsion_points_from_division_polynomial(&curve, 3)
                .expect("exact odd torsion should compute");
            let enumerated_points: Vec<_> = curve
                .points()
                .into_iter()
                .filter(|point| !curve.is_identity(point) && curve.point_order(point) == Some(3))
                .collect();

            prop_assert!(same_point_set::<F17>(&polynomial_points, &enumerated_points));
        }

        #[test]
        fn property_exact_even_torsion_matches_enumeration(
            curve in non_singular_short_weierstrass_curve::<23>(),
        ) {
            let polynomial_points = exact_n_torsion_points_from_division_polynomial(&curve, 6)
                .expect("exact even torsion should compute");
            let enumerated_points: Vec<_> = curve
                .points()
                .into_iter()
                .filter(|point| !curve.is_identity(point) && curve.point_order(point) == Some(6))
                .collect();

            prop_assert!(same_point_set::<F23>(&polynomial_points, &enumerated_points));
        }
    }
}
