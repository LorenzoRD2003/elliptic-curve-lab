use std::collections::BTreeMap;

use crate::fields::{EnumerableFiniteField, SqrtField};

use crate::elliptic_curves::CurveError;

use crate::elliptic_curves::traits::{GroupCurveModel, LiftXCoordinate, PointIndexSampler};

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

fn gcd_usize(mut left: usize, mut right: usize) -> usize {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }

    left
}

fn lcm_usize(left: usize, right: usize) -> usize {
    if left == 0 || right == 0 {
        return 0;
    }

    left / gcd_usize(left, right) * right
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
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem>,
    Self::Point: PartialEq,
{
    /// Returns all finite non-identity points under direct enumeration.
    ///
    /// The current algorithm enumerates every `x` in the base field, lifts the
    /// corresponding points, and deduplicates the `y = 0` case.
    fn finite_points(&self) -> Vec<Self::Point> {
        let mut points = Vec::new();

        for x in Self::BaseField::elements() {
            if let Some((left, right)) = self.points_from_x(x) {
                points.push(left);
                if points.last().is_some_and(|last| *last != right) {
                    points.push(right);
                }
            }
        }

        points
    }

    /// Returns the full point set, with the identity listed first.
    fn points(&self) -> Vec<Self::Point> {
        let mut points = Vec::with_capacity(self.finite_points().len() + 1);
        points.push(self.identity());
        points.extend(self.finite_points());
        points
    }

    /// Returns the total number of points under direct enumeration.
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
    T::BaseField: EnumerableFiniteField<Elem = T::Elem> + SqrtField<Elem = T::Elem>,
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
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem>,
    Self::Point: Clone + PartialEq,
{
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
    T::BaseField: EnumerableFiniteField<Elem = T::Elem> + SqrtField<Elem = T::Elem>,
    T::Point: Clone + PartialEq,
{
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::{
        elliptic_curves::{EnumerableCurveModel, FiniteGroupCurveModel},
        proptest_support::non_singular_short_weierstrass_curve,
    };

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn property_small_enumerable_curves_pass_the_exhaustive_group_axiom_check(
            curve in non_singular_short_weierstrass_curve::<17>(),
        ) {
            prop_assert_eq!(curve.check_group_axioms(), Ok(()));
        }

        #[test]
        fn property_group_structure_matches_order_and_exponent(
            curve in non_singular_short_weierstrass_curve::<17>(),
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
}
