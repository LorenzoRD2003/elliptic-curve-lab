use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::{
    FrobeniusTrace, HasseInterval, HasseMultipleSearchReport, HasseMultipleSearchStep,
    find_annihilating_multiple_in_interval_bsgs, hasse_multiple_search_report,
};
use crate::elliptic_curves::order_from_multiple::mul_scalar_biguint;
use crate::elliptic_curves::traits::{EnumerableCurveModel, GroupCurveModel};
use crate::fields::{EnumerableFiniteField, Field, FiniteField, FiniteFieldDescriptor, SqrtField};
use num_bigint::BigUint;
use std::hash::Hash;

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

/// Curve models with an additive group law that can search one given Hasse interval naively.
///
/// This crate-private extension trait hosts the purely group-theoretic search
/// that starts from one explicit discrete interval `H = [L, U]` and tests
/// `[L]P, [(L+1)]P, ..., [U]P` by one initial scalar setup plus repeated
/// additions of `P`.
pub(crate) trait HasseMultipleSearchCurveModel: GroupCurveModel
where
    Self::Point: Clone,
{
    /// Searches one already-chosen interval from left to right until
    /// `[M]P = O` is found or the interval is exhausted.
    ///
    /// Complexity: one `BigUint` scalar multiplication to build `[L]P`, then
    /// `Θ(|H|)` group additions, where `|H|` is the number of integer
    /// candidates in the supplied interval.
    fn find_annihilating_multiple_in_interval_naive(
        &self,
        point: &Self::Point,
        interval: HasseInterval,
    ) -> Result<HasseMultipleSearchReport<Self::Point>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let lower = interval.lower();
        let upper = interval.upper();
        let mut current = mul_scalar_biguint(self, point, &BigUint::from(lower))?;
        let mut steps = Vec::with_capacity(interval.candidate_count() as usize);
        let mut found = None;

        for candidate_multiple in lower..=upper {
            if candidate_multiple > lower {
                current = self.add(&current, point)?;
            }

            steps.push(HasseMultipleSearchStep::new(
                candidate_multiple,
                current.clone(),
            ));

            if self.is_identity(&current) {
                found = Some(candidate_multiple);
                break;
            }
        }

        Ok(hasse_multiple_search_report(
            interval.q(),
            interval,
            found,
            steps,
        ))
    }

    /// Searches one already-chosen interval with the baby-step/giant-step
    /// method from Algorithm 7.9 in
    /// https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_notes7/
    ///
    /// This helper returns one `M ∈ H(q)` with `[M]P = O`, if found.
    ///
    /// Complexity: Let `c = |H(q) ∩ Z|`. The current implementation chooses
    /// `r = ceil(√c)` and `s = ceil(c/r)`, then performs:
    ///
    /// - `Θ(r)` group additions to build the baby steps
    /// - `Θ(1)` big-scalar multiplications to build `[a]P` and `[r]P`
    /// - `Θ(s)` hash lookups and giant-step additions
    ///
    /// Thus the dominant group-operation count is `Θ(r + s) = Θ(√c)`,
    /// which for Hasse intervals is `Θ(∜q)`.
    fn find_annihilating_multiple_in_interval_bsgs(
        &self,
        point: &Self::Point,
        interval: HasseInterval,
    ) -> Result<Option<u128>, CurveError>
    where
        Self::Point: Eq + Hash,
    {
        find_annihilating_multiple_in_interval_bsgs(self, point, interval)
    }
}

impl<T: GroupCurveModel + ?Sized> HasseMultipleSearchCurveModel for T where T::Point: Clone {}

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
        let base_field = FiniteFieldDescriptor::new(
            Self::BaseField::characteristic(),
            Self::BaseField::extension_degree(),
        )
        .map_err(|_| CurveError::InvalidFrobeniusBaseField {
            characteristic: Self::BaseField::characteristic(),
            extension_degree: Self::BaseField::extension_degree().get(),
        })?;
        let curve_order = self.order() as u64;
        FrobeniusTrace::from_order(base_field, curve_order)
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
    use crate::elliptic_curves::frobenius::relative_frobenius_point;
    use crate::elliptic_curves::{
        AffineCurveModel, RelativeFrobeniusCurveModel, ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};

    type F43 = Fp<43>;

    #[test]
    fn short_weierstrass_relative_frobenius_trait_matches_the_existing_helper() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let point = curve
            .point(F43::zero(), F43::one())
            .expect("(0, 1) should lie on the curve");

        let from_trait = curve
            .relative_frobenius(&point)
            .expect("trait method should evaluate");
        let from_helper =
            relative_frobenius_point(&curve, &point).expect("existing helper should evaluate");

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
}
