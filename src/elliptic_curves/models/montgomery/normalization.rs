use core::fmt;

use crate::elliptic_curves::{AffinePoint, CurveError, MontgomeryCurve};
use crate::fields::traits::{Field, SqrtField};

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

/// Explicit same-field normalization witness from `B y^2 = x^3 + A x^2 + x`
/// to `v^2 = x^3 + A x^2 + x` through `v = √B y`.
pub(crate) struct MontgomeryNormalization<F: Field> {
    source: MontgomeryCurve<F>,
    target: NormalizedMontgomeryCurve<F>,
    sqrt_b: F::Elem,
}

impl<F: Field> MontgomeryNormalization<F> {
    /// Returns the source Montgomery model.
    pub(crate) fn source(&self) -> &MontgomeryCurve<F> {
        &self.source
    }

    /// Returns the normalized target model.
    pub(crate) fn target(&self) -> &NormalizedMontgomeryCurve<F> {
        &self.target
    }

    /// Returns the chosen same-field square root of `B`.
    pub(crate) fn sqrt_b(&self) -> &F::Elem {
        &self.sqrt_b
    }

    /// Transports one source affine point to the normalized target through
    ///
    /// `x -> x`, `v = sqrt(B) y`.
    pub(crate) fn map_source_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.source.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let image = AffinePoint::new(x.clone(), F::mul(&self.sqrt_b, y));
                if !self.target.contains_affine_point(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }
                Ok(image)
            }
        }
    }

    /// Transports one normalized target affine point back to the source model
    /// through `x -> x`, `y = v/√B`.
    pub(crate) fn map_target_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.target.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let original_y = F::div(y, &self.sqrt_b)
                    .expect("sqrt(B) is nonzero on a validated Montgomery curve");
                let image = AffinePoint::new(x.clone(), original_y);
                if !self.source.contains_affine_point(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }
                Ok(image)
            }
        }
    }
}

impl<F: Field> Clone for MontgomeryNormalization<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self {
            source: self.source.clone(),
            target: self.target.clone(),
            sqrt_b: self.sqrt_b.clone(),
        }
    }
}

impl<F: Field + SqrtField> MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    /// Returns the explicit same-field normalization witness
    /// `v = √B y` when `B` is a square in the current base field.
    pub(crate) fn try_normalize(
        &self,
    ) -> Result<MontgomeryNormalization<F>, MontgomeryNormalizationError> {
        let sqrt_b =
            F::sqrt(self.b()).ok_or(MontgomeryNormalizationError::NoSameFieldNormalization)?;
        let target = NormalizedMontgomeryCurve::new(self.a().clone())
            .expect("validated Montgomery curve should normalize to a non-singular B = 1 model");

        Ok(MontgomeryNormalization {
            source: self.clone(),
            target,
            sqrt_b,
        })
    }

    /// Returns the normalized Montgomery companion
    ///
    /// `v^2 = x^3 + A x^2 + x`
    ///
    /// when the same-field scaling witness `sqrt(B)` exists.
    pub fn try_as_normalized_montgomery(
        &self,
    ) -> Result<NormalizedMontgomeryCurve<F>, MontgomeryNormalizationError> {
        self.try_normalize()
            .map(|normalization| normalization.target)
    }
}
