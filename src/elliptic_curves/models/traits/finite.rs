use std::collections::BTreeMap;

use crate::elliptic_curves::{
    CurveError,
    traits::{GroupCurveModel, LiftXCoordinate, LiftedPoints, PointIndexSampler},
};
use crate::fields::traits::EnumerableFiniteField;
use crate::numerics::{lcm_usize, quotients_by_distinct_prime_factors};

/// Structural summary of a small finite abelian curve group.
///
/// For an elliptic curve over a finite field, the group of rational points has
/// the form
///
/// `E(F_q) ≅ Z/n1Z × Z/n2Z`
///
/// with `n1 | n2`. The total order is `n1 * n2`, and the exponent is `n2`.
///
/// This type packages the small-group invariants that the current educational
/// APIs can compute by exhaustive enumeration:
///
/// - `order = #E(F_q)`
/// - `exponent = exp(E(F_q))`
/// - `cyclic`, which is equivalent to the existence of a point of order
///   `#E(F_q)`
/// - `invariant_factors`, which stores `(n1, n2)` only in the non-cyclic case
///
/// When `cyclic` is `true`, the group is isomorphic to `Z/order Z`, so the
/// two-factor decomposition would be the redundant pair `(1, order)`. To keep
/// the surface small and pedagogically direct, this struct uses
/// `invariant_factors = None` in that case.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FiniteAbelianGroupStructure {
    pub order: usize,
    pub exponent: usize,
    pub cyclic: bool,
    pub invariant_factors: Option<(usize, usize)>,
}

fn group_axiom_violation(axiom: &'static str) -> CurveError {
    CurveError::GroupAxiomViolation { axiom }
}

/// Curve models that can be exhaustively enumerated over small finite fields.
///
/// This trait is intentionally narrower than [`super::CurveModel`]. It is
/// meant for educational scenarios where:
///
/// - the base field is small enough to enumerate directly
/// - the curve can be reconstructed by scanning every `x` and lifting
///   `y`-coordinates
///
/// It should not be read as a promise that every curve model in the crate
/// ought to support exhaustive point materialization.
pub trait EnumerableCurveModel: LiftXCoordinate
where
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem>,
    Self::Point: PartialEq,
{
    /// Returns all finite non-identity points under direct enumeration.
    ///
    /// The current algorithm enumerates every `x` in the base field, lifts the
    /// corresponding points, deduplicates the `y = 0` case, and filters out
    /// whichever point the model designates as its identity. That last step is
    /// important for affine models whose neutral element is itself a finite
    /// affine point rather than a distinguished point at infinity.
    fn finite_points(&self) -> Vec<Self::Point> {
        let mut points = Vec::new();

        let mut push_if_nonidentity = |point: Self::Point| {
            if !self.is_identity(&point) {
                points.push(point);
            }
        };

        for x in Self::BaseField::elements() {
            match self
                .lift_x(x)
                .expect("EnumerableCurveModel requires x-fiber lifting to succeed across the whole base field")
            {
                LiftedPoints::NoPoint => {}
                LiftedPoints::OnePoint(point) => push_if_nonidentity(point),
                LiftedPoints::TwoPoints(left, right) => {
                    push_if_nonidentity(left.clone());
                    if left != right {
                        push_if_nonidentity(right);
                    }
                }
            }
        }

        points
    }

    /// Returns the full point set, with the identity listed first.
    ///
    /// This does not assume the identity lies outside the affine chart. Models
    /// with a finite affine neutral element still receive exactly one copy of
    /// that identity here: it is inserted first and filtered out from
    /// [`Self::finite_points`].
    fn points(&self) -> Vec<Self::Point> {
        let mut points = Vec::with_capacity(self.finite_points().len() + 1);
        points.push(self.identity());
        points.extend(self.finite_points());
        points
    }

    /// Returns the total number of points under direct enumeration.
    ///
    /// This is intentionally the exhaustive path: it materializes the full
    /// point set and counts it directly. If another module derives
    /// `#E(F_{q^n})` from Frobenius data, that should stay documented as a
    /// distinct algorithmic route rather than becoming a silent implementation
    /// detail of this method.
    fn order(&self) -> usize {
        self.points().len()
    }

    /// Chooses one point using a minimal index-sampling interface.
    ///
    /// This samples from the fully enumerated point set. It is therefore meant
    /// only for the same small educational settings as [`EnumerableCurveModel`]
    /// itself.
    fn random_point<R>(&self, rng: &mut R) -> Option<Self::Point>
    where
        R: PointIndexSampler,
    {
        let mut points = self.points();
        let index = rng.sample_index(points.len())?;
        if index >= points.len() {
            return None;
        }

        Some(points.swap_remove(index))
    }
}

impl<T> EnumerableCurveModel for T
where
    T: LiftXCoordinate,
    T::BaseField: EnumerableFiniteField<Elem = T::Elem>,
    T::Point: PartialEq,
{
}

/// Small finite curve groups whose full point set can be enumerated honestly.
///
/// This capability is narrower than both [`GroupCurveModel`] and
/// [`EnumerableCurveModel`]. It is meant for educational settings where the
/// whole curve group is small enough that group order and point order can be
/// computed directly.
///
/// In particular, helpers such as [`Self::point_order`] should be read as
/// small-group educational utilities backed by direct traversal, not as
/// efficient generic algorithms for large elliptic-curve groups.
pub trait FiniteGroupCurveModel: GroupCurveModel
where
    Self: EnumerableCurveModel,
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem>,
    Self::Point: Clone + PartialEq,
{
    /// Returns whether `point` has exact order `n`.
    ///
    /// The logic follows the classical exact-order criterion:
    ///
    /// - first check `[n]P = O`
    /// - then rule out every quotient `n / p` where `p` is a prime divisor of `n`
    ///
    /// This is enough because a point has exact order `n` exactly when it is
    /// killed by `[n]` but by no `n / p` for prime `p | n`.
    fn point_has_exact_order(&self, point: &Self::Point, n: usize) -> Result<bool, CurveError> {
        if n == 0 {
            return Err(CurveError::InvalidTorsionOrder { order: n });
        }

        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        if self.is_identity(point) {
            return Ok(n == 1);
        }

        let n_multiple = self.mul_scalar(point, n as u64)?;
        if !self.is_identity(&n_multiple) {
            return Ok(false);
        }

        for quotient in quotients_by_distinct_prime_factors(n) {
            let divisor_multiple = self.mul_scalar(point, quotient as u64)?;
            if self.is_identity(&divisor_multiple) {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Returns the order of a point when the ambient curve group is small
    /// enough to enumerate directly.
    ///
    /// The order is the least positive integer `n` such that `[n]P = O`.
    /// Invalid off-curve inputs return `None` honestly. The identity has order
    /// `1`.
    ///
    /// The current implementation searches by repeated addition up to the full
    /// enumerated group order. That is intentional for small educational
    /// backends and should not be read as an efficient algorithm for large
    /// groups.
    fn point_order(&self, point: &Self::Point) -> Option<usize> {
        if !self.contains(point) {
            return None;
        }

        if self.is_identity(point) {
            return Some(1);
        }

        let mut accumulator = self.identity();
        for n in 1..=self.order() {
            accumulator = self.add(&accumulator, point).ok()?;
            if self.is_identity(&accumulator) {
                return Some(n);
            }
        }

        None
    }

    /// Returns every enumerated point together with its exact order.
    ///
    /// This method is intentionally educational and exhaustive: it computes
    /// the order of each enumerated point by direct traversal.
    fn point_orders(&self) -> Vec<(Self::Point, usize)> {
        self.points()
            .into_iter()
            .map(|point| {
                let order = self
                    .point_order(&point)
                    .expect("enumerated points should have a well-defined order");
                (point, order)
            })
            .collect()
    }

    /// Returns the points whose exact order is `n`.
    ///
    /// The current implementation filters the fully enumerated point set and
    /// is therefore meant only for the same small finite settings as
    /// [`Self::point_order`].
    fn points_of_order(&self, n: usize) -> Vec<Self::Point> {
        self.point_orders()
            .into_iter()
            .filter_map(|(point, order)| (order == n).then_some(point))
            .collect()
    }

    /// Returns the non-identity rational points of exact order `n`.
    ///
    /// This helper filters the fully enumerated rational point set through
    /// [`Self::point_has_exact_order`], so it is honest for both prime and
    /// composite orders in the current small finite educational setting.
    fn points_of_exact_order(&self, n: usize) -> Result<Vec<Self::Point>, CurveError> {
        if n == 0 {
            return Err(CurveError::InvalidTorsionOrder { order: n });
        }

        let mut points = Vec::new();
        for point in self.points() {
            if self.is_identity(&point) {
                continue;
            }

            if self.point_has_exact_order(&point, n)? {
                points.push(point);
            }
        }

        Ok(points)
    }

    /// Returns how many points have each exact order.
    ///
    /// The returned map is keyed by the exact point order and sorted in
    /// ascending order for readable educational output.
    fn order_distribution(&self) -> BTreeMap<usize, usize> {
        let mut distribution = BTreeMap::new();
        for (_, order) in self.point_orders() {
            *distribution.entry(order).or_insert(0) += 1;
        }

        distribution
    }

    /// Returns the exponent of the finite curve group.
    ///
    /// This is the least positive integer `n` such that `[n]P = O` for every
    /// point `P` on the curve. For a finite abelian group, this is the least
    /// common multiple of all point orders.
    fn exponent(&self) -> usize {
        self.order_distribution().into_keys().fold(1, lcm_usize)
    }

    /// Returns one generator when the curve group is cyclic.
    ///
    /// The current implementation searches the enumerated point set for a
    /// point whose order equals the full group order.
    fn generator(&self) -> Option<Self::Point> {
        let group_order = self.order();
        self.point_orders()
            .into_iter()
            .find_map(|(point, order)| (order == group_order).then_some(point))
    }

    /// Returns whether the enumerated finite curve group is cyclic.
    fn is_cyclic(&self) -> bool {
        self.generator().is_some()
    }

    /// Exhaustively checks the promised finite-group laws on the enumerated curve.
    ///
    /// This helper is intentionally educational and brute-force. It is meant
    /// for the same small settings as [`Self::point_order`] and verifies the
    /// current curve model against the following facts:
    ///
    /// - every enumerated point lies on the curve
    /// - the identity point appears in the enumeration
    /// - `P + O = O + P = P`
    /// - `P + (-P) = (-P) + P = O`
    /// - `P + Q = Q + P`
    /// - `(P + Q) + R = P + (Q + R)`
    /// - `[n + m]P = [n]P + [m]P`
    /// - `[ord(P)]P = O`
    /// - `ord(P)` divides `#E(F_q)`
    ///
    /// For the scalar law, the implementation only needs to test
    /// `0 <= n, m < ord(P)`. In a finite group, scalar multiplication by `P`
    /// depends only on the residue class modulo `ord(P)`, so that range is
    /// already exhaustive for the chosen point.
    ///
    /// This is a confidence-building validator for tiny curves, not a
    /// replacement for mathematical proofs or efficient large-group tests.
    fn check_group_axioms(&self) -> Result<(), CurveError> {
        let points = self.points();
        let identity = self.identity();
        let group_order = points.len();

        if !points.iter().all(|point| self.contains(point)) {
            return Err(group_axiom_violation(
                "every enumerated point should lie on the curve",
            ));
        }

        if !points.iter().any(|point| self.is_identity(point)) {
            return Err(group_axiom_violation(
                "the identity point should appear in the enumeration",
            ));
        }

        for point in &points {
            let left_identity_sum = self.add(point, &identity)?;
            if left_identity_sum != *point {
                return Err(group_axiom_violation("P + O = P"));
            }

            let right_identity_sum = self.add(&identity, point)?;
            if right_identity_sum != *point {
                return Err(group_axiom_violation("O + P = P"));
            }

            let inverse = self.neg(point);
            if !self.contains(&inverse) {
                return Err(group_axiom_violation(
                    "the additive inverse of an enumerated point should stay on the curve",
                ));
            }

            if self.add(point, &inverse)? != identity {
                return Err(group_axiom_violation("P + (-P) = O"));
            }

            if self.add(&inverse, point)? != identity {
                return Err(group_axiom_violation("(-P) + P = O"));
            }

            let point_order = self.point_order(point).ok_or_else(|| {
                group_axiom_violation("every enumerated point should have a well-defined order")
            })?;

            if !group_order.is_multiple_of(point_order) {
                return Err(group_axiom_violation("ord(P) divides #E(F_q)"));
            }

            let point_order_u64 = u64::try_from(point_order)
                .expect("tiny educational group orders should fit in u64");

            if self.mul_scalar(point, point_order_u64)? != identity {
                return Err(group_axiom_violation("[ord(P)]P = O"));
            }

            for n in 0..point_order_u64 {
                for m in 0..point_order_u64 {
                    let sum_multiple = self.mul_scalar(point, n + m)?;
                    let split_multiple =
                        self.add(&self.mul_scalar(point, n)?, &self.mul_scalar(point, m)?)?;

                    if sum_multiple != split_multiple {
                        return Err(group_axiom_violation("[n + m]P = [n]P + [m]P"));
                    }
                }
            }
        }

        for left in &points {
            for right in &points {
                let left_plus_right = self.add(left, right)?;
                if !self.contains(&left_plus_right) {
                    return Err(group_axiom_violation(
                        "the sum of enumerated points should stay on the curve",
                    ));
                }

                let right_plus_left = self.add(right, left)?;
                if left_plus_right != right_plus_left {
                    return Err(group_axiom_violation("P + Q = Q + P"));
                }

                for third in &points {
                    let left_grouped = self.add(&left_plus_right, third)?;
                    let right_grouped = self.add(left, &self.add(right, third)?)?;

                    if left_grouped != right_grouped {
                        return Err(group_axiom_violation("(P + Q) + R = P + (Q + R)"));
                    }
                }
            }
        }

        Ok(())
    }

    /// Returns a structural summary of the finite abelian curve group.
    ///
    /// For an elliptic curve over a finite field, the group of rational points
    /// is finite abelian and has the form
    ///
    /// `E(F_q) ≅ Z/n1Z × Z/n2Z`.
    ///
    /// If `N = #E(F_q)` and `m = exp(E(F_q))`, then:
    ///
    /// - when `m = N`, the group is cyclic and isomorphic to `Z/NZ`
    /// - otherwise, one must have `N = n1 * m` with `n1 | m`, so the group is
    ///   `Z/n1Z × Z/mZ`
    ///
    /// The current implementation computes `N` by exhaustive enumeration and
    /// `m` as the least common multiple of all point orders, then reconstructs
    /// the non-cyclic invariant factors as `(N / m, m)`. The divisibility
    /// check `N / m | m` is asserted because it should always hold for finite
    /// elliptic-curve groups.
    fn group_structure(&self) -> FiniteAbelianGroupStructure {
        let order = self.order();
        let exponent = self.exponent();
        let cyclic = exponent == order;

        if order == 1 {
            return FiniteAbelianGroupStructure {
                order,
                exponent,
                cyclic: true,
                invariant_factors: None,
            };
        }

        assert!(
            exponent > 0,
            "finite non-trivial curve groups must have positive exponent"
        );

        if cyclic {
            return FiniteAbelianGroupStructure {
                order,
                exponent,
                cyclic,
                invariant_factors: None,
            };
        }

        assert!(
            order.is_multiple_of(exponent),
            "elliptic-curve group order should be divisible by its exponent"
        );
        let first_factor = order / exponent;
        assert!(
            exponent.is_multiple_of(first_factor),
            "finite elliptic-curve groups should satisfy n1 | n2 in their invariant factors"
        );

        FiniteAbelianGroupStructure {
            order,
            exponent,
            cyclic,
            invariant_factors: Some((first_factor, exponent)),
        }
    }

    /// Returns a small educational description of the finite curve group.
    ///
    /// Finite elliptic-curve groups over finite fields are abelian and admit a
    /// decomposition with at most two cyclic factors. This helper reports that
    /// invariant-factor shape when it can be read directly from the group
    /// order and exponent.
    fn describe_group_structure(&self) -> String {
        let structure = self.group_structure();

        if structure.order == 1 {
            return "trivial group".to_string();
        }

        if structure.cyclic {
            return format!("Z/{}Z", structure.order);
        }

        match structure.invariant_factors {
            Some((left, right)) => format!("Z/{left}Z x Z/{right}Z"),
            None => format!("order {}, exponent {}", structure.order, structure.exponent),
        }
    }
}

impl<T> FiniteGroupCurveModel for T
where
    T: GroupCurveModel + EnumerableCurveModel,
    T::BaseField: EnumerableFiniteField<Elem = T::Elem>,
    T::Point: Clone + PartialEq,
{
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::elliptic_curves::affine::AffinePoint;
    use crate::elliptic_curves::error::CurveError;
    use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
    use crate::elliptic_curves::traits::{
        AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
        LiftXCoordinate, LiftedPoints,
    };
    use crate::fields::{Fp, traits::Field};
    use crate::proptest_support::config::CurveStrategyConfig;
    use crate::proptest_support::elliptic_curves::arb_nonsingular_curve;

    type F3 = Fp<3>;
    type F5 = Fp<5>;
    type F7 = Fp<7>;
    type F41 = Fp<41>;
    type FiniteAffineIdentityPoint = AffinePoint<F3>;
    type Curve5 = ShortWeierstrassCurve<F5>;
    type Curve7 = ShortWeierstrassCurve<F7>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    #[derive(Clone, Debug)]
    struct FiniteAffineIdentityCurve;

    impl FiniteAffineIdentityCurve {
        fn identity_point(&self) -> FiniteAffineIdentityPoint {
            AffinePoint::new(F3::zero(), F3::one())
        }

        fn generator(&self) -> FiniteAffineIdentityPoint {
            AffinePoint::new(F3::one(), F3::from_i64(-1))
        }
    }

    impl CurveModel for FiniteAffineIdentityCurve {
        type Elem = <F3 as Field>::Elem;
        type BaseField = F3;
        type Point = FiniteAffineIdentityPoint;

        fn identity(&self) -> Self::Point {
            self.identity_point()
        }

        fn is_identity(&self, point: &Self::Point) -> bool {
            point == &self.identity_point()
        }

        fn contains(&self, point: &Self::Point) -> bool {
            point == &self.identity_point() || point == &self.generator()
        }
    }

    impl AffineCurveModel for FiniteAffineIdentityCurve {
        fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
            let point = AffinePoint::new(x, y);
            if self.contains(&point) {
                Ok(point)
            } else {
                Err(CurveError::PointNotOnCurve)
            }
        }
    }

    impl LiftXCoordinate for FiniteAffineIdentityCurve {
        fn lift_x(&self, x: Self::Elem) -> Result<LiftedPoints<Self::Point>, CurveError> {
            if F3::eq(&x, &F3::zero()) {
                return Ok(LiftedPoints::OnePoint(self.identity_point()));
            }
            if F3::eq(&x, &F3::one()) {
                return Ok(LiftedPoints::OnePoint(self.generator()));
            }

            Ok(LiftedPoints::NoPoint)
        }
    }

    impl GroupCurveModel for FiniteAffineIdentityCurve {
        fn neg(&self, point: &Self::Point) -> Self::Point {
            if self.is_identity(point) {
                self.identity()
            } else {
                self.generator()
            }
        }

        fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
            if !self.contains(left) || !self.contains(right) {
                return Err(CurveError::PointNotOnCurve);
            }
            if self.is_identity(left) {
                return Ok(right.clone());
            }
            if self.is_identity(right) {
                return Ok(left.clone());
            }

            Ok(self.identity())
        }
    }

    fn f5_noncyclic_curve() -> Curve5 {
        Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
    }

    fn f7_curve() -> Curve7 {
        Curve7::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn finite_points_excludes_a_finite_affine_identity() {
        let curve = FiniteAffineIdentityCurve;

        assert_eq!(curve.finite_points(), vec![curve.generator()]);
    }

    #[test]
    fn points_includes_a_finite_affine_identity_exactly_once() {
        let curve = FiniteAffineIdentityCurve;

        assert_eq!(
            curve.points(),
            vec![curve.identity_point(), curve.generator()]
        );
    }

    #[test]
    fn shared_small_group_helpers_respect_models_with_finite_affine_identity() {
        let curve = FiniteAffineIdentityCurve;

        assert_eq!(curve.order(), 2);
        assert_eq!(curve.point_order(&curve.identity()), Some(1));
        assert_eq!(curve.point_order(&curve.generator()), Some(2));
        assert_eq!(curve.points_of_order(1), vec![curve.identity()]);
        assert_eq!(curve.points_of_exact_order(2), Ok(vec![curve.generator()]));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn property_small_enumerable_curves_pass_the_exhaustive_group_axiom_check(
            curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default()),
        ) {
            prop_assert_eq!(curve.check_group_axioms(), Ok(()));
        }

        #[test]
        fn property_group_structure_matches_order_and_exponent(
            curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default()),
        ) {
            let structure = curve.group_structure();

            prop_assert_eq!(structure.order, curve.order());
            prop_assert_eq!(structure.exponent, curve.exponent());
            if let Some((left, right)) = structure.invariant_factors {
                prop_assert_eq!(left * right, structure.order);
                prop_assert!(right.is_multiple_of(left));
            } else {
                prop_assert!(structure.cyclic);
            }
        }
    }

    #[test]
    fn point_has_exact_order_rejects_zero_degree() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(
            curve.point_has_exact_order(&point, 0),
            Err(CurveError::InvalidTorsionOrder { order: 0 })
        );
    }

    #[test]
    fn point_has_exact_order_rejects_points_outside_the_curve() {
        let curve = f41_curve();
        let off_curve = AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

        assert_eq!(
            curve.point_has_exact_order(&off_curve, 2),
            Err(CurveError::PointNotOnCurve)
        );
    }

    #[test]
    fn point_has_exact_order_accepts_exact_order_two_point() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(curve.point_has_exact_order(&point, 2), Ok(true));
    }

    #[test]
    fn point_has_exact_order_rejects_points_killed_by_a_prime_factor_quotient() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(curve.point_has_exact_order(&point, 4), Ok(false));
    }

    #[test]
    fn points_of_exact_order_rejects_zero_degree() {
        let curve = f41_curve();

        assert_eq!(
            curve.points_of_exact_order(0),
            Err(CurveError::InvalidTorsionOrder { order: 0 })
        );
    }

    #[test]
    fn points_of_exact_order_returns_non_identity_points_of_exact_order() {
        let curve = f41_curve();
        let points = curve
            .points_of_exact_order(2)
            .expect("degree two should be valid");

        assert_eq!(
            points,
            vec![curve.point(F41::from_i64(40), F41::from_i64(0)).unwrap()]
        );
        assert!(points.iter().all(|point| !curve.is_identity(point)));
    }

    #[test]
    fn finds_nonzero_two_torsion_points() {
        let curve = f5_noncyclic_curve();
        let points = curve
            .points_of_exact_order(2)
            .expect("degree two should be valid");

        assert_eq!(
            points,
            vec![
                curve.point(F5::from_i64(0), F5::from_i64(0)).unwrap(),
                curve.point(F5::from_i64(1), F5::from_i64(0)).unwrap(),
                curve.point(F5::from_i64(4), F5::from_i64(0)).unwrap(),
            ]
        );
    }

    #[test]
    fn identity_is_not_reported_as_nontrivial_exact_order_point() {
        let curve = f5_noncyclic_curve();
        let identity = curve.identity();

        assert_eq!(curve.point_has_exact_order(&identity, 2), Ok(false));

        let points = curve
            .points_of_exact_order(2)
            .expect("degree two should be valid");
        assert!(!points.contains(&identity));
    }

    #[test]
    fn points_of_exact_order_matches_known_small_order_three_example() {
        let curve = f7_curve();
        let points = curve
            .points_of_exact_order(3)
            .expect("degree three should be valid");

        assert_eq!(
            points,
            vec![
                curve.point(F7::from_i64(3), F7::from_i64(1)).unwrap(),
                curve.point(F7::from_i64(3), F7::from_i64(6)).unwrap(),
            ]
        );
    }
}
