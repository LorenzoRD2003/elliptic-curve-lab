use crate::elliptic_curves::{
    CurveError,
    frobenius::hasse::search::HasseIntervalSearchCurveModel,
    frobenius::{
        FrobeniusTrace,
        hasse::{HasseBoundReport, HasseMultipleSearchReport},
    },
    traits::{EnumerableCurveModel, GroupCurveModel},
};
use crate::fields::traits::*;
use crate::fields::{
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{EnumerableFiniteField, FiniteField, SqrtField},
};

/// Curve models over a finite field that expose the relative Frobenius `π_q`.
///
/// If the base field has size `q`, this trait models the curve endomorphism
/// induced by the coordinate action
///
/// `π_q(x, y) = (x^q, y^q)`.
///
/// For a curve defined over `F_q`, this relative Frobenius sends points of the
/// model back to the same curve. The current trait is intentionally narrow:
/// it records only the point-level endomorphism needed by later educational
/// helpers such as characteristic-equation checks.
///
/// Complexity:
/// The trait does not prescribe one asymptotic cost for
/// [`RelativeFrobeniusCurveModel::relative_frobenius`]. Concrete
/// implementations should document that cost. The default squared helper below
/// simply doubles the cost of one relative-Frobenius evaluation.
pub trait RelativeFrobeniusCurveModel: GroupCurveModel
where
    Self::BaseField: FiniteField,
    Self::Point: Clone,
{
    /// Applies the relative Frobenius `π_q` to a point on the curve.
    ///
    /// Implementations should return [`CurveError::PointNotOnCurve`] for
    /// invalid off-curve inputs.
    fn relative_frobenius(&self, point: &Self::Point) -> Result<Self::Point, CurveError>;

    /// Applies the square `π_q^2` of the relative Frobenius.
    ///
    /// The default implementation reuses [`Self::relative_frobenius`] twice.
    ///
    /// Complexity:
    /// Two calls to [`Self::relative_frobenius`].
    fn relative_frobenius_squared(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        let first = self.relative_frobenius(point)?;
        self.relative_frobenius(&first)
    }
}

/// Enumerable curve models that can recover the Frobenius trace over `F_q`.
///
/// This capability is intentionally stronger than
/// [`RelativeFrobeniusCurveModel`]: in addition to exposing the point-level
/// relative Frobenius, the model must live in a small finite setting where the
/// full point set can be enumerated honestly and `#E(F_q)` can therefore be
/// recovered by direct counting.
///
/// Complexity:
/// The trace helper below is intentionally exhaustive. Its dominant cost is
/// the curve's full rational-point enumeration.
pub trait FrobeniusTraceCurveModel: EnumerableCurveModel
where
    Self::BaseField:
        EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem> + FiniteField,
    Self::Point: Clone + PartialEq,
{
    /// Computes the Frobenius trace from an exhaustive point count on `E(F_q)`.
    ///
    /// The current implementation is intentionally exhaustive and therefore
    /// intended only for small enumerable finite fields, but the counting
    /// identity itself is the general finite-field formula
    ///
    /// `t = q + 1 - #E(F_q)`.
    ///
    /// Complexity:
    /// In the current implementation this is the cost of enumerating all
    /// rational points of the curve, plus `Θ(1)` integer post-processing.
    fn frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        let characteristic = Self::BaseField::characteristic().to_biguint();
        let base_field =
            FiniteFieldDescriptor::new(characteristic, Self::BaseField::extension_degree())
                .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                    characteristic: Self::BaseField::characteristic().to_biguint(),
                    extension_degree: Self::BaseField::extension_degree().get(),
                })?;
        let curve_order = self.order();
        FrobeniusTrace::from_order(base_field, curve_order)
    }

    /// Verifies Hasse's bound `|t| <= 2 sqrt(q)` through the equivalent exact
    /// integer inequality `t^2 <= 4q`.
    ///
    /// Complexity: the cost of [`Self::frobenius_trace`] plus `Θ(1)` integer
    /// arithmetic.
    fn verify_hasse_bound(&self) -> Result<HasseBoundReport, CurveError> {
        self.frobenius_trace()
            .map(HasseBoundReport::from_frobenius_trace)
    }

    /// Searches the discrete Hasse interval `H(q)` from left to right until
    /// it finds the first `M` with `[M]P = O`.
    ///
    /// The current educational implementation is the naive route from the
    /// notes:
    ///
    /// 1. compute the initial image `[L]P`, where `L = min H(q)`
    /// 2. step through the interval by repeated addition of `P`
    /// 3. stop at the first identity image
    ///
    /// Complexity: One `BigUint` scalar multiplication to build `[L]P`,
    /// followed by `Θ(|H(q)|)` group additions. Since `|H(q)| = Θ(√q)`,
    /// this is a `Θ(√q)`-addition search after the initial setup.
    fn find_annihilating_multiple_in_hasse_interval_naive(
        &self,
        point: &Self::Point,
    ) -> Result<HasseMultipleSearchReport<Self::Point>, CurveError>
    where
        Self: GroupCurveModel,
    {
        let trace = self.frobenius_trace()?;
        self.find_annihilating_multiple_in_interval_naive(point, trace.hasse_interval())
    }
}

impl<T> FrobeniusTraceCurveModel for T
where
    T: EnumerableCurveModel,
    T::BaseField: EnumerableFiniteField<Elem = T::Elem> + SqrtField<Elem = T::Elem> + FiniteField,
    T::Point: Clone + PartialEq,
{
}

#[cfg(test)]
mod tests {

    use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
    use crate::elliptic_curves::traits::{AffineCurveModel, RelativeFrobeniusCurveModel};
    use crate::fields::traits::Field;

    type F43 = crate::fields::Fp43;

    #[test]
    fn short_weierstrass_relative_frobenius_trait_matches_the_existing_helper() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("(0, 1) should lie on the curve");

        let from_trait = curve
            .relative_frobenius(&point)
            .expect("trait method should evaluate");
        let from_helper = curve
            .relative_frobenius_point(&point)
            .expect("existing helper should evaluate");

        assert_eq!(from_trait, from_helper);
    }

    #[test]
    fn default_relative_frobenius_squared_reuses_the_relative_map_twice() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("(0, 1) should lie on the curve");

        let squared = curve
            .relative_frobenius_squared(&point)
            .expect("default square should evaluate");
        let once = curve
            .relative_frobenius(&point)
            .expect("first relative Frobenius should evaluate");
        let twice = curve
            .relative_frobenius(&once)
            .expect("second relative Frobenius should evaluate");

        assert_eq!(squared, twice);
    }

    #[test]
    fn relative_frobenius_keeps_base_field_points_fixed() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("(0, 1) should lie on the curve");

        let image = curve
            .relative_frobenius(&point)
            .expect("relative Frobenius should evaluate");

        assert_eq!(image, point);
    }
}
