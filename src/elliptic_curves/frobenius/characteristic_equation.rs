use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::{FrobeniusCharacteristicPolynomial, FrobeniusTrace};
use crate::elliptic_curves::traits::{FrobeniusTraceCurveModel, RelativeFrobeniusCurveModel};
use crate::fields::{EnumerableFiniteField, Field, FiniteField, FiniteFieldDescriptor, SqrtField};

/// Pointwise verification data for the Frobenius characteristic equation.
///
/// For a curve over `F_q` with Frobenius trace `t`, the relative Frobenius
/// `π_q` satisfies `π_q^2(P) - [t]π_q(P) + [q]P = O`.
///
/// In the current educational finite-field setting, this struct records the
/// concrete pointwise terms of that identity on a point of the represented model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicEquationCheck<P> {
    point: P,
    pi_q: P,
    pi_q_squared: P,
    trace_term: P,
    field_order_term: P,
    lhs: P,
    holds: bool,
}

impl<P> FrobeniusCharacteristicEquationCheck<P> {
    /// Returns the input point `P`.
    pub fn point(&self) -> &P {
        &self.point
    }

    /// Returns `π_q(P)`.
    pub fn pi_q(&self) -> &P {
        &self.pi_q
    }

    /// Returns `π_q^2(P)`.
    pub fn pi_q_squared(&self) -> &P {
        &self.pi_q_squared
    }

    /// Returns `[t]π_q(P)`.
    pub fn trace_term(&self) -> &P {
        &self.trace_term
    }

    /// Returns `[q]P`.
    pub fn field_order_term(&self) -> &P {
        &self.field_order_term
    }

    /// Returns the left-hand side `π_q^2(P) - [t]π_q(P) + [q]P`.
    pub fn lhs(&self) -> &P {
        &self.lhs
    }

    /// Returns whether the characteristic equation holds at the checked point.
    pub fn holds(&self) -> bool {
        self.holds
    }
}

/// Exhaustive verification report for the Frobenius characteristic equation.
///
/// This report is intended for the same tiny finite-field educational setting
/// as [`FrobeniusTraceCurveModel`]. It records the Frobenius trace used for
/// the check, how many rational points were tested, and the full pointwise
/// checks that failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicEquationExhaustiveReport<P> {
    frobenius_trace: FrobeniusTrace,
    checked_points: usize,
    failed_checks: Vec<FrobeniusCharacteristicEquationCheck<P>>,
}

impl<P> FrobeniusCharacteristicEquationExhaustiveReport<P> {
    /// Returns the Frobenius trace package used by the exhaustive check.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the number of checked rational points.
    pub fn checked_points(&self) -> usize {
        self.checked_points
    }

    /// Returns the pointwise checks that failed.
    pub fn failed_checks(&self) -> &[FrobeniusCharacteristicEquationCheck<P>] {
        &self.failed_checks
    }

    /// Returns whether the characteristic equation held on every checked point.
    pub fn all_hold(&self) -> bool {
        self.failed_checks.is_empty()
    }
}

/// Verifies the characteristic equation of the relative Frobenius at one point.
///
/// If `χ_{π_q}(T) = T^2 - tT + q` is the characteristic polynomial of the
/// relative Frobenius on a curve over `F_q`, this helper computes the pointwise
/// terms of `π_q^2(P) - [t]π_q(P) + [q]P`
/// and reports whether the result is the identity.
///
/// Complexity:
/// This helper performs:
/// - two relative-Frobenius evaluations
/// - one signed scalar multiplication by `t`
/// - one scalar multiplication by `q`
/// - one subtraction and one addition in the curve group
///
/// With the current `GroupCurveModel` defaults, the dominant cost
/// is the two scalar multiplications, so the group-law complexity is
/// `Θ(log |t| + log q)` additions/doublings, plus the cost of the two
/// relative-Frobenius evaluations and one membership check.
pub fn verify_frobenius_characteristic_equation_at_point<E: RelativeFrobeniusCurveModel>(
    curve: &E,
    point: &E::Point,
    characteristic_polynomial: &FrobeniusCharacteristicPolynomial,
) -> Result<FrobeniusCharacteristicEquationCheck<E::Point>, CurveError>
where
    E::BaseField: FiniteField,
    E::Point: Clone + PartialEq,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let curve_base_field = curve_base_field::<E>();
    let polynomial_base_field = characteristic_polynomial.base_field();
    if &curve_base_field != polynomial_base_field {
        return Err(CurveError::IncompatibleFrobeniusBaseField {
            curve_characteristic: curve_base_field.characteristic,
            curve_extension_degree: curve_base_field.extension_degree.get(),
            polynomial_characteristic: polynomial_base_field.characteristic,
            polynomial_extension_degree: polynomial_base_field.extension_degree.get(),
        });
    }

    let pi_q = curve.relative_frobenius(point)?;
    let pi_q_squared = curve.relative_frobenius_squared(point)?;
    let trace_term = curve.mul_scalar_signed(&pi_q, characteristic_polynomial.trace())?;
    let q_scalar = u64::try_from(characteristic_polynomial.field_order()).map_err(|_| {
        CurveError::UnsupportedFrobeniusFieldOrder {
            field_order: characteristic_polynomial.field_order(),
        }
    })?;
    let field_order_term = curve.mul_scalar(point, q_scalar)?;
    let lhs_without_q = curve.sub(&pi_q_squared, &trace_term)?;
    let lhs = curve.add(&lhs_without_q, &field_order_term)?;
    let holds = curve.is_identity(&lhs);

    Ok(FrobeniusCharacteristicEquationCheck {
        point: point.clone(),
        pi_q,
        pi_q_squared,
        trace_term,
        field_order_term,
        lhs,
        holds,
    })
}

/// Verifies the Frobenius characteristic equation on every rational point of a
/// small enumerable curve over `F_q`.
///
/// The current implementation first computes the Frobenius trace from the
/// exhaustive point count, derives `χ_{π_q}(T) = T^2 - tT + q`, and then
/// reuses [`verify_frobenius_characteristic_equation_at_point`] on every point
/// returned by the curve's enumeration surface.
///
/// Complexity: If `N = #E(F_q)`, this helper performs:
/// - one exhaustive Frobenius-trace recovery by counting
/// - one full enumeration of the rational point set
/// - `N` pointwise characteristic-equation checks
///
/// So its total cost is the cost of exhaustive point counting and enumeration,
/// plus `N` times the pointwise verification cost above. In the current small
/// educational setting, this should be read as an explicitly exhaustive
/// `Θ(N (log |t| + log q))` group-law pass on top of the enumeration work.
pub fn verify_frobenius_characteristic_equation_exhaustive<E>(
    curve: &E,
) -> Result<FrobeniusCharacteristicEquationExhaustiveReport<E::Point>, CurveError>
where
    E: FrobeniusTraceCurveModel + RelativeFrobeniusCurveModel,
    E::BaseField: EnumerableFiniteField<Elem = E::Elem> + SqrtField<Elem = E::Elem> + FiniteField,
    E::Point: Clone + PartialEq,
{
    let frobenius_trace = curve.frobenius_trace()?;
    let characteristic_polynomial = frobenius_trace.characteristic_polynomial();
    let points = curve.points();
    let checked_points = points.len();
    let mut failed_checks = Vec::new();

    for point in points {
        let check = verify_frobenius_characteristic_equation_at_point(
            curve,
            &point,
            &characteristic_polynomial,
        )?;
        if !check.holds() {
            failed_checks.push(check);
        }
    }

    Ok(FrobeniusCharacteristicEquationExhaustiveReport {
        frobenius_trace,
        checked_points,
        failed_checks,
    })
}

fn curve_base_field<E: RelativeFrobeniusCurveModel>() -> FiniteFieldDescriptor
where
    E::BaseField: FiniteField,
{
    FiniteFieldDescriptor::new(
        E::BaseField::characteristic(),
        E::BaseField::extension_degree(),
    )
    .expect("finite field implementations should expose internally consistent metadata")
}
