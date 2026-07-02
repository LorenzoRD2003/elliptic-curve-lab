use crate::elliptic_curves::traits::{
    CurveModel, FiniteGroupCurveModel, GroupCurveModel, ScalarInput,
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    error::{IsogenyConstructionError, IsogenyError},
    kernel::{KernelDescription, ReducedKernelDescription},
    traits::Isogeny,
};
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};

/// Scalar-multiplication isogeny `[n] : E -> E` on a small explicit curve group.
///
/// For a non-zero integer `n`, elliptic-curve multiplication by `n`
///
/// `[n](P) = P + P + ... + P`
///
/// is an example of an isogeny from a curve to itself.
///
/// In this educational implementation:
///
/// - the domain and codomain are the same curve value
/// - the degree is reported as `n^2`
/// - `kernel_description()` exposes the currently honest kernel description
/// - in reduced small-field cases, `kernel_points()` still recovers the
///   visible rational points killed by `[n]`
pub struct ScalarMultiplicationIsogeny<C: GroupCurveModel> {
    curve: C,
    scalar: BigUint,
    kernel_points: Vec<C::Point>,
}

impl<C: FiniteGroupCurveModel> ScalarMultiplicationIsogeny<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    /// Builds the scalar-multiplication isogeny `[n]`.
    ///
    /// The current constructor is intentionally restricted to small finite
    /// curve groups so it can materialize the kernel
    ///
    /// `E[n] = { P in E(F_q) : [n]P = O }`.
    ///
    /// Scalar `0` is rejected because this crate reserves the isogeny surface
    /// for the usual non-constant multiplication-by-`n` maps.
    pub fn new(curve: C, scalar: impl ScalarInput) -> Result<Self, IsogenyError> {
        let scalar = scalar.into_biguint_scalar();
        if scalar.is_zero() {
            return Err(IsogenyConstructionError::ZeroScalarIsNotIsogeny.into());
        }

        let identity = curve.identity();
        let kernel_points = curve
            .points()
            .into_iter()
            .map(|point| -> Result<_, IsogenyError> {
                let image = curve.mul_scalar(&point, &scalar)?;
                Ok((point, image == identity))
            })
            .collect::<Result<Vec<_>, IsogenyError>>()?
            .into_iter()
            .filter_map(|(point, kills_point)| kills_point.then_some(point))
            .collect();

        Ok(Self {
            curve,
            scalar,
            kernel_points,
        })
    }

    /// Returns the underlying scalar `n` in `[n]`.
    pub fn scalar(&self) -> &BigUint {
        &self.scalar
    }

    pub fn curve(&self) -> &C {
        &self.curve
    }

    pub fn kernel_points(&self) -> &[<C as CurveModel>::Point] {
        &self.kernel_points
    }
}

impl<C: FiniteGroupCurveModel> Isogeny<C, C> for ScalarMultiplicationIsogeny<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    fn domain(&self) -> &C {
        &self.curve
    }

    fn codomain(&self) -> &C {
        &self.curve
    }

    fn degree(&self) -> usize {
        (&self.scalar * &self.scalar)
            .to_usize()
            .expect("educational scalar-multiplication degrees should fit in usize")
    }

    fn evaluate(&self, p: &C::Point) -> Result<C::Point, IsogenyError> {
        self.curve.mul_scalar(p, &self.scalar).map_err(Into::into)
    }

    fn kernel_description(&self) -> KernelDescription<C> {
        let factorization = self.scalar_characteristic_factorization();
        if factorization.p_power_exponent() == 0 {
            KernelDescription::Reduced(
                ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints {
                    points: self.kernel_points.clone(),
                    degree: self.kernel_points.len(),
                },
            )
        } else {
            KernelDescription::Mixed(self.mixed_kernel_description().expect(
                "visible reduced kernel enumeration should succeed on stored finite curves",
            ))
        }
    }

    fn kernel_points(&self) -> Vec<C::Point> {
        self.kernel_points.clone()
    }
}
