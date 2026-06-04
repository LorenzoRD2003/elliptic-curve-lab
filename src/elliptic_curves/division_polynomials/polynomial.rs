use std::collections::BTreeMap;

use crate::{
    elliptic_curves::{AffinePoint, CurveError, ShortWeierstrassCurve, traits::CurveModel},
    fields::Field,
    polynomials::{DensePolynomial, evaluation::evaluate_dense},
};

use super::DivisionPolynomialError;

/// For a short-Weierstrass curve `E: y^2 = x^3 + ax + b`,
/// the first division polynomials are:
///
/// - `ψ_0 = 0`
/// - `ψ_1 = 1`
/// - `ψ_2 = 2y`
/// - `ψ_3 = 3x^4 + 6ax^2 + 12bx - a^2`
/// - `ψ_4 = 4y (x^6 + 5ax^4 + 20bx^3 - 5a^2x^2 - 4abx - 8b^2 - a^3)`
///
/// The important structural distinction is:
/// - for odd `n`, `ψ_n ∈ F[x]`
/// - for even `n`, `ψ_n ∈ yF[x]`
///
/// This enum therefore records whether the object lives directly in `F[x]`
/// or in `yF[x]`, while reusing `DensePolynomial<F>` for the polynomial part
/// in the `x`-coordinate.
#[derive(Clone, Debug)]
pub enum DivisionPolynomialForm<F: Field> {
    /// A division polynomial that lies in `F[x]`.
    InX(DensePolynomial<F>),
    /// A division polynomial of the form `y * f(x)` with `f(x) ∈ F[x]`.
    YTimes(DensePolynomial<F>),
}

/// Backward-compatible alias for the division-polynomial shape enum.
pub type DivisionPolynomial<F> = DivisionPolynomialForm<F>;

/// The `x`-coordinate criterion used when checking whether a division
/// polynomial vanishes through `x` alone.
///
/// For odd indices, the full division polynomial already lies in `F[x]`, so
/// the relevant criterion is just `ψ_n(x)`.
///
/// For even indices, `ψ_n = y ε_n(x)` is not itself an `x`-polynomial, so the
/// meaningful `x`-only criterion is the stripped factor `ε_n(x)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DivisionPolynomialXCriterionKind {
    /// Use the odd division polynomial `ψ_n(x)`.
    OddDivisionPolynomial,
    /// Use the even `x`-factor `ε_n(x)` in `ψ_n = y ε_n(x)`.
    EvenYStrippedFactor,
}

impl<F: Field> PartialEq for DivisionPolynomialForm<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::InX(lhs), Self::InX(rhs)) => lhs == rhs,
            (Self::YTimes(lhs), Self::YTimes(rhs)) => lhs == rhs,
            _ => false,
        }
    }
}

impl<F: Field> DivisionPolynomialForm<F> {
    /// Wraps a polynomial that lies directly in `F[x]`.
    pub fn x_polynomial(polynomial: DensePolynomial<F>) -> Self {
        Self::InX(polynomial)
    }

    /// Wraps a polynomial factor `f(x)` to represent `y * f(x)`.
    pub fn y_times_x_polynomial(polynomial: DensePolynomial<F>) -> Self {
        Self::YTimes(polynomial)
    }

    /// Returns whether this division polynomial lies directly in `F[x]`.
    pub fn is_x_polynomial(&self) -> bool {
        matches!(self, Self::InX(_))
    }

    /// Returns whether this division polynomial has the form `y * f(x)`.
    pub fn is_y_times_x_polynomial(&self) -> bool {
        matches!(self, Self::YTimes(_))
    }

    /// Returns the stored polynomial factor in `F[x]`.
    ///
    /// For `InX(f)`, this returns `f`.
    ///
    /// For `YTimes(f)`, this returns the factor `f` in the
    /// decomposition `y * f(x)`.
    pub fn x_factor(&self) -> &DensePolynomial<F> {
        match self {
            Self::InX(polynomial) | Self::YTimes(polynomial) => polynomial,
        }
    }
}

/// Returns the `x`-coordinate criterion kind attached to index `n`.
///
/// - odd `n` use the honest `x`-polynomial `ψ_n(x)`
/// - even `n` use only the stripped factor `ε_n(x)` from `ψ_n = y ε_n(x)`
///
/// Index `0` is rejected because division polynomials are indexed from `1`
/// upward.
pub fn division_polynomial_x_criterion_kind(
    n: usize,
) -> Result<DivisionPolynomialXCriterionKind, DivisionPolynomialError> {
    if n == 0 {
        Err(DivisionPolynomialError::ZeroIndex)
    } else if n.is_multiple_of(2) {
        Ok(DivisionPolynomialXCriterionKind::EvenYStrippedFactor)
    } else {
        Ok(DivisionPolynomialXCriterionKind::OddDivisionPolynomial)
    }
}

/// Returns the base short-Weierstrass division polynomials `ψ_0` through
/// `ψ_4`.
///
/// This helper currently covers only the explicit low-degree formulas:
///
/// - `ψ_0 = 0`
/// - `ψ_1 = 1`
/// - `ψ_2 = 2y`
/// - `ψ_3 = 3x^4 + 6ax^2 + 12bx - a^2`
/// - `ψ_4 = 4y (x^6 + 5ax^4 + 20bx^3 - 5a^2x^2 - 4abx - 8b^2 - a^3)`
///
/// For larger indices, the current milestone-7 scaffold reports
/// [`DivisionPolynomialError::UnsupportedIndex`].
///
/// Complexity:
///
/// this helper performs only a constant amount of field arithmetic and builds
/// coefficient vectors of bounded size, so its running time is `O(1)`.
pub fn division_polynomial_base<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<DivisionPolynomialForm<F>, DivisionPolynomialError> {
    fn scaled<F: Field>(coefficient: &F::Elem, scalar: i64) -> F::Elem {
        F::mul(coefficient, &F::from_i64(scalar))
    }

    let a = curve.a();
    let b = curve.b();
    let a2 = F::square(a);
    let a3 = F::mul(&a2, a);
    let b2 = F::square(b);
    let ab = F::mul(a, b);

    match n {
        0 => Ok(DivisionPolynomialForm::x_polynomial(DensePolynomial::new(
            Vec::new(),
        ))),
        1 => Ok(DivisionPolynomialForm::x_polynomial(
            DensePolynomial::constant(F::one()),
        )),
        2 => Ok(DivisionPolynomialForm::y_times_x_polynomial(
            DensePolynomial::constant(F::from_i64(2)),
        )),
        3 => Ok(DivisionPolynomialForm::x_polynomial(DensePolynomial::new(
            vec![
                F::neg(&a2),
                scaled::<F>(b, 12),
                scaled::<F>(a, 6),
                F::zero(),
                F::from_i64(3),
            ],
        ))),
        4 => {
            let constant_term = F::add(&scaled::<F>(&b2, -32), &scaled::<F>(&a3, -4));

            Ok(DivisionPolynomialForm::y_times_x_polynomial(
                DensePolynomial::new(vec![
                    constant_term,
                    scaled::<F>(&ab, -16),
                    scaled::<F>(&a2, -20),
                    scaled::<F>(b, 80),
                    scaled::<F>(a, 20),
                    F::zero(),
                    F::from_i64(4),
                ]),
            ))
        }
        _ => Err(DivisionPolynomialError::UnsupportedIndex { n }),
    }
}

/// Returns the short-Weierstrass right-hand side
/// `R(x) = x^3 + ax + b`.
///
/// This polynomial is the algebraic replacement for `y^2` used throughout the
/// recursive division-polynomial formulas below. Whenever a recurrence would
/// naturally produce an even power of `y`, the implementation substitutes
/// `y^2 = R(x)` so the result stays inside `F[x]`.
fn short_weierstrass_rhs_polynomial<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
) -> DensePolynomial<F> {
    DensePolynomial::new(vec![
        curve.b().clone(),
        curve.a().clone(),
        F::zero(),
        F::one(),
    ])
}

/// Squares a dense polynomial using the crate's current naive multiplication.
///
/// With the present dense backend, multiplying degree-`d` polynomials costs
/// `O(d^2)` field operations.
fn dense_square<F: Field>(polynomial: &DensePolynomial<F>) -> DensePolynomial<F> {
    polynomial.mul(polynomial)
}

/// Cubes a dense polynomial by one squaring and one extra multiplication.
///
/// With the current dense backend, this remains `O(d^2)` up to a larger
/// constant factor, since both dense multiplications are quadratic in the
/// degree of the input.
fn dense_cube<F: Field>(polynomial: &DensePolynomial<F>) -> DensePolynomial<F> {
    dense_square(polynomial).mul(polynomial)
}

/// Recursively computes the odd division polynomial `ψ_n ∈ F[x]`.
///
/// - odd indices are represented directly as polynomials in `F[x]`
/// - even indices are represented separately through their `F[x]` factor
///   inside `ψ_n = y ε_n(x)`
///
/// Let `R(x) = x^3 + ax + b`. The implementation uses the classical
/// short-Weierstrass recurrences, but rewrites every even power of `y` through
/// `y^2 = R(x)`.
///
/// Terminology:
/// - `δ_n(x) = ψ_n(x)` for odd `n`
/// - `ε_n(x)` for the unique polynomial with `ψ_n = y ε_n(x)` when `n` is even
///
/// the code computes `δ_{2m+1}` by:
/// - if `m` is even:
///   `δ_{2m+1} = R(x)^2 ε_{m+2} ε_m^3 - δ_{m-1} δ_{m+1}^3`
/// - if `m` is odd:
///   `δ_{2m+1} = δ_{m+2} δ_m^3 - R(x)^2 ε_{m-1} ε_{m+1}^3`
///
/// Memoization: previously computed odd and even subproblems are stored in separate
/// `BTreeMap`s, so each index is expanded at most once during a top-level
/// call. This avoids the exponential recomputation that a direct recursive
/// tree would otherwise incur.
///
/// Complexity: the degree of `ψ_n` is `Θ(n^2)`, so at index `k` the dominant dense
/// multiplications cost `Θ(k^4)` with the current naive backend. Since
/// memoization computes each relevant index only once up to `n`, the total
/// running time for a top-level request is roughly
/// `Σ_{k≤n} Θ(k^4) = Θ(n^5)`, ignoring lower-order cache lookup costs.
///
/// Future work: with FFT, the same recursive and memoized schema would have
/// computational complexity of Θ(n^3 log n).
fn odd_division_polynomial_inner<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    odd_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
    even_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
    if let Some(polynomial) = odd_cache.get(&n) {
        return Ok(polynomial.clone());
    }

    let polynomial = match n {
        0 => return Err(DivisionPolynomialError::ZeroIndex),
        1 | 3 => match division_polynomial_base(curve, n)? {
            DivisionPolynomial::InX(polynomial) => polynomial,
            DivisionPolynomial::YTimes(_) => {
                unreachable!("odd base case must lie in F[x]")
            }
        },
        _ if n.is_multiple_of(2) => {
            return Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n });
        }
        _ => {
            let m = (n - 1) / 2;
            let rhs_squared = dense_square(&short_weierstrass_rhs_polynomial(curve));

            if m.is_multiple_of(2) {
                let even_m_plus_2 =
                    even_division_polynomial_inner(curve, m + 2, odd_cache, even_cache)?;
                let even_m = even_division_polynomial_inner(curve, m, odd_cache, even_cache)?;
                let odd_m_minus_1 =
                    odd_division_polynomial_inner(curve, m - 1, odd_cache, even_cache)?;
                let odd_m_plus_1 =
                    odd_division_polynomial_inner(curve, m + 1, odd_cache, even_cache)?;

                rhs_squared
                    .mul(&even_m_plus_2.mul(&dense_cube(&even_m)))
                    .sub(&odd_m_minus_1.mul(&dense_cube(&odd_m_plus_1)))
            } else {
                let odd_m_plus_2 =
                    odd_division_polynomial_inner(curve, m + 2, odd_cache, even_cache)?;
                let odd_m = odd_division_polynomial_inner(curve, m, odd_cache, even_cache)?;
                let even_m_minus_1 =
                    even_division_polynomial_inner(curve, m - 1, odd_cache, even_cache)?;
                let even_m_plus_1 =
                    even_division_polynomial_inner(curve, m + 1, odd_cache, even_cache)?;

                odd_m_plus_2
                    .mul(&dense_cube(&odd_m))
                    .sub(&rhs_squared.mul(&even_m_minus_1.mul(&dense_cube(&even_m_plus_1))))
            }
        }
    };

    odd_cache.insert(n, polynomial.clone());
    Ok(polynomial)
}

/// Recursively computes the `F[x]` factor `ε_n(x)` for an even division
/// polynomial `ψ_n = y ε_n(x)`.
///
/// Mathematical normalization:
///
/// instead of trying to represent even division polynomials as elements of
/// `F[x]` directly, the implementation strips the unavoidable `y` factor and
/// stores only `ε_n(x)`. This is the mathematically honest reason why the
/// public `even_division_polynomial` helper returns the factor in `F[x]`
/// rather than a fictitious polynomial for `ψ_n` itself.
///
/// Using the same notation as above, the classical even recurrence
/// is rewritten into parity-aware `F[x]` formulas:
///
/// - if `m` is even:
///   `ε_{2m} = (1/2) ε_m (ε_{m+2} δ_{m-1}^2 - ε_{m-2} δ_{m+1}^2)`
/// - if `m` is odd:
///   `ε_{2m} = (1/2) δ_m (δ_{m+2} ε_{m-1}^2 - δ_{m-2} ε_{m+1}^2)`
///
/// The factor `1/2` is valid because validated short-Weierstrass curves in
/// this crate exclude characteristic `2`, so the element `2` is invertible.
///
/// Memoization and complexity:
///
/// as with odd indices, subresults are cached by index in `BTreeMap`s. The
/// asymptotic cost is again dominated by dense polynomial multiplication and
/// is therefore `Θ(n^5)` for a top-level request under the current naive
/// multiplication backend.
fn even_division_polynomial_inner<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    odd_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
    even_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
    if let Some(polynomial) = even_cache.get(&n) {
        return Ok(polynomial.clone());
    }

    let polynomial = match n {
        0 => return Err(DivisionPolynomialError::ZeroIndex),
        2 | 4 => match division_polynomial_base(curve, n)? {
            DivisionPolynomial::YTimes(polynomial) => polynomial,
            DivisionPolynomial::InX(_) => unreachable!("even base case must lie in yF[x]"),
        },
        _ if n % 2 == 1 => return Err(DivisionPolynomialError::UnsupportedIndex { n }),
        _ => {
            let m = n / 2;
            let inverse_two = F::inverse(&F::from_i64(2))
                .expect("validated short-Weierstrass curves have characteristic different from 2");

            if m.is_multiple_of(2) {
                let even_m = even_division_polynomial_inner(curve, m, odd_cache, even_cache)?;
                let even_m_plus_2 =
                    even_division_polynomial_inner(curve, m + 2, odd_cache, even_cache)?;
                let even_m_minus_2 =
                    even_division_polynomial_inner(curve, m - 2, odd_cache, even_cache)?;
                let odd_m_minus_1 =
                    odd_division_polynomial_inner(curve, m - 1, odd_cache, even_cache)?;
                let odd_m_plus_1 =
                    odd_division_polynomial_inner(curve, m + 1, odd_cache, even_cache)?;

                even_m
                    .mul(
                        &even_m_plus_2
                            .mul(&dense_square(&odd_m_minus_1))
                            .sub(&even_m_minus_2.mul(&dense_square(&odd_m_plus_1))),
                    )
                    .scale(&inverse_two)
            } else {
                let odd_m = odd_division_polynomial_inner(curve, m, odd_cache, even_cache)?;
                let odd_m_plus_2 =
                    odd_division_polynomial_inner(curve, m + 2, odd_cache, even_cache)?;
                let odd_m_minus_2 =
                    odd_division_polynomial_inner(curve, m - 2, odd_cache, even_cache)?;
                let even_m_minus_1 =
                    even_division_polynomial_inner(curve, m - 1, odd_cache, even_cache)?;
                let even_m_plus_1 =
                    even_division_polynomial_inner(curve, m + 1, odd_cache, even_cache)?;

                odd_m
                    .mul(
                        &odd_m_plus_2
                            .mul(&dense_square(&even_m_minus_1))
                            .sub(&odd_m_minus_2.mul(&dense_square(&even_m_plus_1))),
                    )
                    .scale(&inverse_two)
            }
        }
    };

    even_cache.insert(n, polynomial.clone());
    Ok(polynomial)
}

/// Returns the odd-index division polynomial `ψ_n` as an element of `F[x]`.
///
/// Precondition:
///
/// - `n` must be odd
/// - `n >= 1`
///
/// If `n` is even, this helper reports
/// [`DivisionPolynomialError::EvenIndexRequiresYFactor`] because even
/// division polynomials live in `yF[x]`, not in `F[x]` itself.
///
/// Algorithm:
///
/// this function dispatches to the recursive parity-aware recurrence described
/// on [`odd_division_polynomial_inner`], while memoizing odd and even
/// subproblems in separate caches for the duration of the call.
///
/// Complexity:
///
/// with the current dense naive multiplication backend, the total cost is
/// roughly `Θ(n^5)` for index `n`. If the crate later gains an FFT-based
/// multiplication routine so that multiplying degree-`d` polynomials costs
/// `O(d log d)` instead of `O(d^2)`, the same memoized recurrence strategy
/// would improve to about `Θ(n^3 log n)`.
pub fn odd_division_polynomial<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
    let mut odd_cache = BTreeMap::new();
    let mut even_cache = BTreeMap::new();
    odd_division_polynomial_inner(curve, n, &mut odd_cache, &mut even_cache)
}

/// Returns the polynomial factor `f_n(x) ∈ F[x]` for an even-index division
/// polynomial `ψ_n = y · f_n(x)`.
///
/// Precondition:
///
/// - `n` must be even
/// - `n >= 2`
///
/// Educational note:
///
/// this helper returns only the `F[x]` factor. The full even division
/// polynomial is obtained by multiplying the result by `y`.
///
/// Algorithm:
///
/// this function computes the normalized factor `ε_n(x)` through the
/// parity-aware recursive formulas documented on
/// [`even_division_polynomial_inner`], again with per-call memoization of both
/// odd and even indices.
///
/// Complexity:
///
/// with the current dense naive multiplication backend, the total cost is
/// roughly `Θ(n^5)` for index `n`. With FFT-based polynomial multiplication,
/// the same high-level recursion and memoization strategy would drop to about
/// `Θ(n^3 log n)`.
pub fn even_division_polynomial<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
    let mut odd_cache = BTreeMap::new();
    let mut even_cache = BTreeMap::new();
    even_division_polynomial_inner(curve, n, &mut odd_cache, &mut even_cache)
}

/// Returns the short-Weierstrass division polynomial in its honest public
/// shape.
///
/// The current dense polynomial type is parameterized by the field family `F`
/// rather than by the raw coefficient type `F::Elem`, so this API returns
/// `DivisionPolynomialForm<F>`. Semantically, this is the intended public
/// surface:
///
/// - odd `n` map to `DivisionPolynomialForm::InX(...)`
/// - even `n` map to `DivisionPolynomialForm::YTimes(...)`
pub fn division_polynomial<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<DivisionPolynomialForm<F>, DivisionPolynomialError> {
    match n {
        0..=4 => division_polynomial_base(curve, n),
        _ if n.is_multiple_of(2) => Ok(DivisionPolynomialForm::YTimes(even_division_polynomial(
            curve, n,
        )?)),
        _ => Ok(DivisionPolynomialForm::InX(odd_division_polynomial(
            curve, n,
        )?)),
    }
}

/// Evaluates an odd division polynomial `ψ_n ∈ F[x]` at an `x`-coordinate.
/// - `n` must be odd
/// - `n >= 1`
pub fn evaluate_odd_division_polynomial_at_x<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    x: &F::Elem,
) -> Result<F::Elem, DivisionPolynomialError> {
    let polynomial = odd_division_polynomial(curve, n)?;
    evaluate_dense(&polynomial, x).map_err(Into::into)
}

/// Evaluates the `F[x]` factor `ε_n(x)` in an even division polynomial
/// `ψ_n = y ε_n(x)`.
/// - `n` must be even
/// - `n >= 2`
pub fn evaluate_even_division_polynomial_factor_at_x<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    x: &F::Elem,
) -> Result<F::Elem, DivisionPolynomialError> {
    let polynomial = even_division_polynomial(curve, n)?;
    evaluate_dense(&polynomial, x).map_err(Into::into)
}

/// Evaluates the division-polynomial `x`-criterion attached to index `n`.
///
/// This is the right helper when callers only have an `x`-coordinate:
///
/// - if `n` is odd, it evaluates `ψ_n(x)`
/// - if `n` is even, it evaluates the stripped factor `ε_n(x)` from
///   `ψ_n = y ε_n(x)`
pub fn evaluate_division_polynomial_x_criterion<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    x: &F::Elem,
) -> Result<F::Elem, DivisionPolynomialError> {
    match division_polynomial_x_criterion_kind(n)? {
        DivisionPolynomialXCriterionKind::OddDivisionPolynomial => {
            evaluate_odd_division_polynomial_at_x(curve, n, x)
        }
        DivisionPolynomialXCriterionKind::EvenYStrippedFactor => {
            evaluate_even_division_polynomial_factor_at_x(curve, n, x)
        }
    }
}

/// Evaluates `ψ_n(P)` at a finite affine point `P`.
///
/// - if `n` is odd, the implementation evaluates `ψ_n(x)` at the point's
///   `x`-coordinate
/// - if `n` is even, it evaluates the factor `ε_n(x)` and then multiplies by
///   the point's `y`-coordinate
///
/// The current implementation does not support evaluating at the point at infinity.
pub fn evaluate_division_polynomial_at_point<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
    point: &AffinePoint<F>,
) -> Result<F::Elem, DivisionPolynomialError> {
    if !curve.contains(point) {
        return Err(DivisionPolynomialError::Curve(CurveError::PointNotOnCurve));
    }

    let AffinePoint::Finite { x, y } = point else {
        return Err(DivisionPolynomialError::PointAtInfinityNotSupported);
    };

    let criterion = evaluate_division_polynomial_x_criterion(curve, n, x)?;

    match division_polynomial_x_criterion_kind(n)? {
        DivisionPolynomialXCriterionKind::OddDivisionPolynomial => Ok(criterion),
        DivisionPolynomialXCriterionKind::EvenYStrippedFactor => Ok(F::mul(y, &criterion)),
    }
}

#[cfg(test)]
mod tests {
    use super::{
        DivisionPolynomial, DivisionPolynomialForm, DivisionPolynomialXCriterionKind,
        division_polynomial, division_polynomial_base, division_polynomial_x_criterion_kind,
        evaluate_division_polynomial_at_point, evaluate_division_polynomial_x_criterion,
        evaluate_even_division_polynomial_factor_at_x, evaluate_odd_division_polynomial_at_x,
        even_division_polynomial, odd_division_polynomial,
    };
    use crate::{
        elliptic_curves::{AffineCurveModel, AffinePoint, GroupCurveModel, ShortWeierstrassCurve},
        fields::{Field, Fp},
        polynomials::DensePolynomial,
    };

    type F17 = Fp<17>;
    type F23 = Fp<23>;

    #[test]
    fn x_polynomial_variant_reports_itself_honestly() {
        let polynomial = DensePolynomial::<F17>::constant(F17::one());
        let division_polynomial = DivisionPolynomial::x_polynomial(polynomial.clone());

        assert!(division_polynomial.is_x_polynomial());
        assert!(!division_polynomial.is_y_times_x_polynomial());
        assert_eq!(division_polynomial.x_factor(), &polynomial);
    }

    #[test]
    fn y_times_x_polynomial_variant_reports_itself_honestly() {
        let polynomial = DensePolynomial::<F17>::new(vec![F17::elem_from_u64(2), F17::one()]);
        let division_polynomial = DivisionPolynomial::y_times_x_polynomial(polynomial.clone());

        assert!(division_polynomial.is_y_times_x_polynomial());
        assert!(!division_polynomial.is_x_polynomial());
        assert_eq!(division_polynomial.x_factor(), &polynomial);
    }

    #[test]
    fn division_polynomial_base_covers_psi_zero_through_four() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
            .expect("curve should be non-singular");

        let psi_0 = division_polynomial_base(&curve, 0).expect("psi_0 should exist");
        let psi_1 = division_polynomial_base(&curve, 1).expect("psi_1 should exist");
        let psi_2 = division_polynomial_base(&curve, 2).expect("psi_2 should exist");
        let psi_3 = division_polynomial_base(&curve, 3).expect("psi_3 should exist");
        let psi_4 = division_polynomial_base(&curve, 4).expect("psi_4 should exist");

        assert_eq!(
            psi_0,
            DivisionPolynomial::x_polynomial(DensePolynomial::<F17>::new(Vec::new()))
        );
        assert_eq!(
            psi_1,
            DivisionPolynomial::x_polynomial(DensePolynomial::<F17>::constant(F17::one()))
        );
        assert_eq!(
            psi_2,
            DivisionPolynomial::y_times_x_polynomial(DensePolynomial::<F17>::constant(
                F17::from_i64(2)
            ))
        );
        assert_eq!(
            psi_3,
            DivisionPolynomial::x_polynomial(DensePolynomial::<F17>::new(vec![
                F17::from_i64(-4),
                F17::from_i64(36),
                F17::from_i64(12),
                F17::zero(),
                F17::from_i64(3),
            ]))
        );
        assert_eq!(
            psi_4,
            DivisionPolynomial::y_times_x_polynomial(DensePolynomial::<F17>::new(vec![
                F17::from_i64(-320),
                F17::from_i64(-96),
                F17::from_i64(-80),
                F17::from_i64(240),
                F17::from_i64(40),
                F17::zero(),
                F17::from_i64(4),
            ]))
        );
    }

    #[test]
    fn odd_division_polynomial_matches_base_cases_and_rejects_even_indices() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            odd_division_polynomial(&curve, 1).expect("psi_1 should exist"),
            DensePolynomial::<F17>::constant(F17::one())
        );
        assert_eq!(
            odd_division_polynomial(&curve, 3).expect("psi_3 should exist"),
            DensePolynomial::<F17>::new(vec![
                F17::from_i64(-4),
                F17::from_i64(36),
                F17::from_i64(12),
                F17::zero(),
                F17::from_i64(3),
            ])
        );
        assert_eq!(
            odd_division_polynomial(&curve, 0),
            Err(super::DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            odd_division_polynomial(&curve, 2),
            Err(super::DivisionPolynomialError::EvenIndexRequiresYFactor { n: 2 })
        );
    }

    #[test]
    fn even_division_polynomial_matches_base_cases_and_rejects_odd_indices() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            even_division_polynomial(&curve, 2).expect("psi_2 factor should exist"),
            DensePolynomial::<F17>::constant(F17::from_i64(2))
        );
        assert_eq!(
            even_division_polynomial(&curve, 4).expect("psi_4 factor should exist"),
            DensePolynomial::<F17>::new(vec![
                F17::from_i64(-320),
                F17::from_i64(-96),
                F17::from_i64(-80),
                F17::from_i64(240),
                F17::from_i64(40),
                F17::zero(),
                F17::from_i64(4),
            ])
        );
        assert_eq!(
            even_division_polynomial(&curve, 0),
            Err(super::DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            even_division_polynomial(&curve, 3),
            Err(super::DivisionPolynomialError::UnsupportedIndex { n: 3 })
        );
    }

    #[test]
    fn public_division_polynomial_api_uses_the_expected_shape_by_parity() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
            .expect("curve should be non-singular");

        assert_eq!(
            division_polynomial(&curve, 3).expect("psi_3 should exist"),
            DivisionPolynomialForm::InX(DensePolynomial::<F17>::new(vec![
                F17::from_i64(-4),
                F17::from_i64(36),
                F17::from_i64(12),
                F17::zero(),
                F17::from_i64(3),
            ]))
        );
        assert_eq!(
            division_polynomial(&curve, 4).expect("psi_4 should exist"),
            DivisionPolynomialForm::YTimes(DensePolynomial::<F17>::new(vec![
                F17::from_i64(-320),
                F17::from_i64(-96),
                F17::from_i64(-80),
                F17::from_i64(240),
                F17::from_i64(40),
                F17::zero(),
                F17::from_i64(4),
            ]))
        );
    }

    #[test]
    fn x_criterion_kind_dispatches_by_parity() {
        assert_eq!(
            division_polynomial_x_criterion_kind(0),
            Err(super::DivisionPolynomialError::ZeroIndex)
        );
        assert_eq!(
            division_polynomial_x_criterion_kind(3),
            Ok(DivisionPolynomialXCriterionKind::OddDivisionPolynomial)
        );
        assert_eq!(
            division_polynomial_x_criterion_kind(4),
            Ok(DivisionPolynomialXCriterionKind::EvenYStrippedFactor)
        );
    }

    #[test]
    fn x_criterion_evaluation_matches_the_existing_odd_and_even_helpers() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let x = F23::elem_from_u64(1);

        assert_eq!(
            evaluate_division_polynomial_x_criterion(&curve, 5, &x).unwrap(),
            evaluate_odd_division_polynomial_at_x(&curve, 5, &x).unwrap()
        );
        assert_eq!(
            evaluate_division_polynomial_x_criterion(&curve, 6, &x).unwrap(),
            evaluate_even_division_polynomial_factor_at_x(&curve, 6, &x).unwrap()
        );
    }

    #[test]
    fn recursive_odd_evaluation_at_x_matches_point_evaluation() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let point = AffinePoint::<F23>::new(F23::elem_from_u64(1), F23::elem_from_u64(12));

        let x_value = evaluate_odd_division_polynomial_at_x(&curve, 5, &F23::elem_from_u64(1))
            .expect("psi_5(x) should evaluate");
        let point_value = evaluate_division_polynomial_at_point(&curve, 5, &point)
            .expect("psi_5(P) should evaluate");
        let negated_value = evaluate_division_polynomial_at_point(&curve, 5, &point.neg())
            .expect("psi_5(-P) should evaluate");

        assert!(F23::eq(&x_value, &point_value));
        assert!(F23::eq(&point_value, &negated_value));
    }

    #[test]
    fn recursive_even_factor_evaluation_matches_point_evaluation_up_to_y() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let point = AffinePoint::<F23>::new(F23::elem_from_u64(1), F23::elem_from_u64(12));

        let factor =
            evaluate_even_division_polynomial_factor_at_x(&curve, 6, &F23::elem_from_u64(1))
                .expect("epsilon_6(x) should evaluate");
        let point_value = evaluate_division_polynomial_at_point(&curve, 6, &point)
            .expect("psi_6(P) should evaluate");
        let negated_value = evaluate_division_polynomial_at_point(&curve, 6, &point.neg())
            .expect("psi_6(-P) should evaluate");
        let expected = F23::mul(&F23::elem_from_u64(12), &factor);

        assert!(F23::eq(&expected, &point_value));
        assert!(F23::eq(&F23::neg(&point_value), &negated_value));
    }

    #[test]
    fn evaluate_at_identity_is_handled_explicitly() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let infinity = AffinePoint::<F23>::infinity();

        assert_eq!(
            evaluate_division_polynomial_at_point(&curve, 3, &infinity),
            Err(super::DivisionPolynomialError::PointAtInfinityNotSupported)
        );
    }

    #[test]
    fn point_evaluation_rejects_off_curve_inputs() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let off_curve = AffinePoint::<F23>::new(F23::zero(), F23::zero());

        assert_eq!(
            evaluate_division_polynomial_at_point(&curve, 3, &off_curve),
            Err(super::DivisionPolynomialError::Curve(
                crate::elliptic_curves::CurveError::PointNotOnCurve
            ))
        );
    }

    #[test]
    fn odd_division_polynomial_vanishes_on_known_three_torsion_points() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let positive = curve
            .point(F23::elem_from_u64(8), F23::elem_from_u64(5))
            .expect("sample three-torsion point should lie on the curve");
        let negative = positive.neg();

        assert!(curve.is_torsion_point(&positive, 3));
        assert!(curve.is_torsion_point(&negative, 3));

        let positive_value = evaluate_division_polynomial_at_point(&curve, 3, &positive)
            .expect("psi_3 should evaluate on the sample three-torsion point");
        let negative_value = evaluate_division_polynomial_at_point(&curve, 3, &negative)
            .expect("psi_3 should evaluate on the negated three-torsion point");

        assert!(F23::eq(&positive_value, &F23::zero()));
        assert!(F23::eq(&negative_value, &F23::zero()));
    }

    #[test]
    fn odd_division_polynomial_does_not_vanish_on_generic_non_torsion_point() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::elem_from_u64(2), F23::elem_from_u64(3))
            .expect("curve should be non-singular");
        let point = curve
            .point(F23::elem_from_u64(1), F23::elem_from_u64(12))
            .expect("sample point should lie on the curve");

        assert!(!curve.is_torsion_point(&point, 3));

        let value = evaluate_division_polynomial_at_point(&curve, 3, &point)
            .expect("psi_3 should evaluate on a generic non-torsion point");

        assert!(!F23::eq(&value, &F23::zero()));
    }
}
