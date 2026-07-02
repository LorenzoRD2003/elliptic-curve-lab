use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::{AffinePoint, CurveError, MontgomeryCurve};
use crate::fields::traits::SqrtField;

/// Error returned when the requested same-field Montgomery normalization does
/// not exist.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MontgomeryNormalizationError {
    /// The scaling witness `sqrt(B)` does not exist over the current base field.
    NoSameFieldNormalization,
}

impl fmt::Display for MontgomeryNormalizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoSameFieldNormalization => write!(
                formatter,
                "Montgomery normalization v = sqrt(B) y is not available over the current base field"
            ),
        }
    }
}

impl std::error::Error for MontgomeryNormalizationError {}

/// Normalized Montgomery model  `v^2 = x^3 + A x^2 + x`.
pub struct NormalizedMontgomeryCurve<F: Field> {
    a: F::Elem,
}

impl<F: Field> NormalizedMontgomeryCurve<F> {
    /// Builds one validated normalized Montgomery descriptor.
    pub fn new(a: F::Elem) -> Result<Self, CurveError> {
        MontgomeryCurve::<F>::new(a.clone(), F::one())?;
        Ok(Self { a })
    }

    /// Returns the coefficient `A`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the cubic right-hand side `x^3 + A x^2 + x`.
    pub fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        F::add(&F::add(&F::cube(x), &F::mul(self.a(), &F::square(x))), x)
    }

    /// Returns whether one affine point satisfies the normalized equation.
    pub fn contains_affine_point(&self, point: &AffinePoint<F>) -> bool {
        match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => F::eq(&F::square(y), &self.rhs_value(x)),
        }
    }

    /// Builds one checked affine point on the normalized model.
    pub fn point(&self, x: F::Elem, y: F::Elem) -> Result<AffinePoint<F>, CurveError> {
        let point = AffinePoint::new(x, y);
        if self.contains_affine_point(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }

    /// Returns the normalized model as the ambient Montgomery descriptor with
    /// `B = 1`.
    pub fn as_montgomery_curve(&self) -> MontgomeryCurve<F>
    where
        F::Elem: Clone,
    {
        MontgomeryCurve::new(self.a.clone(), F::one())
            .expect("validated normalized Montgomery model should stay non-singular")
    }

    /// Returns the defining equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("v^2 = x^3 + ({})x^2 + x", self.a)
    }
}

impl<F: Field> Clone for NormalizedMontgomeryCurve<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self { a: self.a.clone() }
    }
}

impl<F: Field> PartialEq for NormalizedMontgomeryCurve<F>
where
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a() == other.a()
    }
}

impl<F: Field> Eq for NormalizedMontgomeryCurve<F> where F::Elem: Eq {}

impl<F: Field> fmt::Display for NormalizedMontgomeryCurve<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for NormalizedMontgomeryCurve<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("NormalizedMontgomeryCurve")
            .field(
                "equation",
                &format_args!("v^2 = x^3 + ({:?})x^2 + x", self.a()),
            )
            .field("a", self.a())
            .finish()
    }
}

impl<F: Field + SqrtField> MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    /// Returns the normalized Montgomery companion
    ///
    /// `v^2 = x^3 + A x^2 + x`
    ///
    /// when the same-field scaling witness `sqrt(B)` exists.
    pub fn try_as_normalized_montgomery(
        &self,
    ) -> Result<NormalizedMontgomeryCurve<F>, MontgomeryNormalizationError> {
        if F::sqrt(self.b()).is_none() {
            return Err(MontgomeryNormalizationError::NoSameFieldNormalization);
        }

        NormalizedMontgomeryCurve::new(self.a().clone())
            .map_err(|_| MontgomeryNormalizationError::NoSameFieldNormalization)
    }
}
