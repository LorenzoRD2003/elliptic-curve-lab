use std::collections::HashMap;
use std::hash::Hash;

use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::frobenius::FrobeniusOrbit;
use crate::elliptic_curves::frobenius::absolute_frobenius_power_point;
use crate::elliptic_curves::frobenius::orbit::{
    orbit_from_successor_by_key, partition_point_orbits_by_key,
};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::torsion::points_of_exact_order;
use crate::elliptic_curves::traits::{FiniteGroupCurveModel, RelativeFrobeniusCurveModel};
use crate::fields::{EnumerableFiniteField, FiniteField, SqrtField};

/// Frobenius data for one rational point of exact order `n`.
///
/// The same point-level value object is used for both:
///
/// - the relative Frobenius `π_q` on `E(F_q)`, where the result is currently
///   tautological
/// - the absolute Frobenius `π_p^k` on torsion points represented in an
///   extension field, where nontrivial motion can already appear
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusOnExactTorsionPoint<P> {
    point: P,
    frobenius_image: P,
    minimal_absolute_frobenius_fixing_power: Option<u32>,
}

impl<P: PartialEq> FrobeniusOnExactTorsionPoint<P> {
    /// Returns the torsion point being checked.
    pub fn point(&self) -> &P {
        &self.point
    }

    /// Returns the Frobenius image of the checked point.
    pub fn frobenius_image(&self) -> &P {
        &self.frobenius_image
    }

    /// Returns whether the applied Frobenius map fixes the point.
    ///
    /// Complexity: `Θ(1)`.
    pub fn fixed_by_frobenius(&self) -> bool {
        self.point == self.frobenius_image
    }

    /// Returns the smallest positive `d` such that `π_p^d(P) = P`, when this
    /// metadata was computed by an absolute-Frobenius helper.
    ///
    /// Relative-Frobenius reports do not populate this invariant and therefore
    /// return `None`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn minimal_absolute_frobenius_fixing_power(&self) -> Option<u32> {
        self.minimal_absolute_frobenius_fixing_power
    }

    /// Returns whether `π_p^k(P) = P`, when absolute-Frobenius fixing metadata
    /// is available.
    ///
    /// Since the smallest fixing power `d` satisfies `π_p^d(P) = P`, the point
    /// is fixed by `π_p^k` exactly when `d | k`. We also treat `k = 0` as the
    /// identity iterate.
    ///
    /// Relative-Frobenius reports do not populate this invariant and therefore
    /// return `None`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn fixed_by_absolute_frobenius_power(&self, power: u32) -> Option<bool> {
        self.minimal_absolute_frobenius_fixing_power
            .map(|minimal_power| {
                if power == 0 {
                    true
                } else {
                    power.is_multiple_of(minimal_power)
                }
            })
    }
}

/// Report for a Frobenius action on the exact-`n` rational torsion.
///
/// The current implementation records the non-identity points of exact order
/// `n` found in the represented finite-field model and their images under the
/// Frobenius map chosen by the calling helper.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusOnExactTorsionReport<P> {
    exact_order: usize,
    points: Vec<FrobeniusOnExactTorsionPoint<P>>,
}

impl<P: PartialEq> FrobeniusOnExactTorsionReport<P> {
    /// Returns the exact torsion order `n`.
    pub fn exact_order(&self) -> usize {
        self.exact_order
    }

    /// Returns the checked rational exact-`n` torsion points and their images.
    pub fn points(&self) -> &[FrobeniusOnExactTorsionPoint<P>] {
        &self.points
    }

    /// Returns whether every listed torsion point is fixed by the applied Frobenius map.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of listed exact-`n` torsion points.
    pub fn all_fixed(&self) -> bool {
        self.points.iter().all(|point| point.fixed_by_frobenius())
    }

    /// Returns how many listed torsion points are fixed by the applied Frobenius map.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of listed exact-`n` torsion points.
    pub fn fixed_count(&self) -> usize {
        self.points
            .iter()
            .filter(|point| point.fixed_by_frobenius())
            .count()
    }

    /// Returns how many listed torsion points are moved by the applied Frobenius map.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of listed exact-`n` torsion points.
    pub fn moved_count(&self) -> usize {
        self.points.len() - self.fixed_count()
    }

    /// Returns how many listed torsion points are already fixed by the prime-field
    /// absolute Frobenius.
    ///
    /// For reports produced by `π_p` on a base-defined curve viewed over
    /// `F_{p^r}`, this is the number of exact-`n` torsion points that already
    /// descend to `F_p`.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of exact-`n` torsion points.
    pub fn prime_field_rational_count(&self) -> usize {
        self.fixed_count()
    }

    /// Returns how many listed torsion points are moved by the prime-field
    /// absolute Frobenius.
    ///
    /// For reports produced by `π_p` on a base-defined curve viewed over
    /// `F_{p^r}`, this is the number of exact-`n` torsion points that are
    /// visible in the chosen extension field but do not yet descend to `F_p`.
    ///
    /// For other Frobenius actions, this is simply an educational alias for
    /// [`Self::moved_count`].
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of exact-`n` torsion points.
    pub fn extension_only_count(&self) -> usize {
        self.moved_count()
    }

    /// Returns how many listed torsion points are fixed by `π_p^k`, when the
    /// report carries absolute-Frobenius fixing metadata.
    ///
    /// Relative-Frobenius reports do not populate this invariant and therefore
    /// return `None`.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of listed exact-`n`
    /// torsion points.
    pub fn fixed_by_absolute_frobenius_power_count(&self, power: u32) -> Option<usize> {
        let mut count = 0;
        for point in &self.points {
            match point.fixed_by_absolute_frobenius_power(power) {
                Some(true) => count += 1,
                Some(false) => {}
                None => return None,
            }
        }
        Some(count)
    }

    /// Returns how many listed torsion points have minimal absolute-Frobenius
    /// fixing power exactly `d`, when that metadata is available.
    ///
    /// Relative-Frobenius reports do not populate this invariant and therefore
    /// return `None`.
    ///
    /// Complexity: `Θ(m)`, where `m` is the number of listed exact-`n`
    /// torsion points.
    pub fn count_with_minimal_absolute_frobenius_fixing_power(&self, power: u32) -> Option<usize> {
        let mut count = 0;
        for point in &self.points {
            match point.minimal_absolute_frobenius_fixing_power() {
                Some(minimal_power) if minimal_power == power => count += 1,
                Some(_) => {}
                None => return None,
            }
        }
        Some(count)
    }

    /// Partitions the reported exact-`n` torsion points into Frobenius orbits
    /// under the same point-to-image action stored in this report.
    ///
    /// This integrates the orbit layer with the torsion layer without asking
    /// the caller to recompute the Frobenius action from the curve again.
    ///
    /// Complexity:
    /// If `m` is the number of listed exact-`n` torsion points, this method
    /// first builds a hashed point-to-image table of size `m` and then
    /// partitions the induced permutation into orbits. Under the usual
    /// constant-time hash-table model, the total bookkeeping is `Θ(m)`.
    pub fn orbits(&self) -> Vec<FrobeniusOrbit<P>>
    where
        P: Clone + Eq + Hash,
    {
        let frobenius_images_by_point: HashMap<P, P> = self
            .points
            .iter()
            .map(|entry| (entry.point().clone(), entry.frobenius_image().clone()))
            .collect();

        partition_point_orbits_by_key(
            self.points.iter().map(|entry| entry.point().clone()),
            |point| point.clone(),
            |point| {
                orbit_from_successor_by_key(
                    point.clone(),
                    self.points.len(),
                    |orbit_point| orbit_point.clone(),
                    |current| {
                        Ok(frobenius_images_by_point
                            .get(current)
                            .expect("orbit report should contain the image source point")
                            .clone())
                    },
                )
            },
        )
        .expect("stored torsion report should induce valid Frobenius orbits")
    }

    /// Returns how many Frobenius orbits appear in the reported exact-`n`
    /// torsion set.
    ///
    /// Complexity: `Θ(m)` under the usual constant-time hash-table model.
    pub fn orbit_count(&self) -> usize
    where
        P: Clone + Eq + Hash,
    {
        self.orbits().len()
    }

    /// Returns the periods of all Frobenius orbits in the reported exact-`n`
    /// torsion set.
    ///
    /// Complexity: `Θ(m)` under the usual constant-time hash-table model.
    pub fn orbit_periods(&self) -> Vec<usize>
    where
        P: Clone + Eq + Hash,
    {
        self.orbits()
            .into_iter()
            .map(|orbit| orbit.period())
            .collect()
    }
}

/// Computes the relative Frobenius action on the rational points of exact
/// order `n`.
///
/// This helper first enumerates the non-identity rational points of exact
/// order `n` and then applies the relative Frobenius `π_q` to each one.
///
/// In the current API surface, those points already lie in `E(F_q)`, so the
/// result is mathematically tautological: every reported point is fixed by
/// `π_q`. Even so, this makes the action explicit and gives later extension
/// work a stable report shape to build on.
///
/// Complexity:
/// If `N = #E(F_q)` and `m = #E(F_q)[n]_{exact}`, this helper performs:
/// - one exhaustive exact-order torsion scan across the rational point set
/// - `m` relative-Frobenius evaluations
///
/// So the total cost is the current cost of [`points_of_exact_order`], plus
/// `Θ(m)` relative-Frobenius evaluations.
pub fn relative_frobenius_on_exact_torsion<E: FiniteGroupCurveModel + RelativeFrobeniusCurveModel>(
    curve: &E,
    exact_order: usize,
) -> Result<FrobeniusOnExactTorsionReport<E::Point>, CurveError>
where
    E::BaseField: EnumerableFiniteField<Elem = E::Elem> + SqrtField<Elem = E::Elem> + FiniteField,
    E::Point: Clone + PartialEq,
{
    let torsion_points = points_of_exact_order(curve, exact_order)?;
    let mut points = Vec::with_capacity(torsion_points.len());

    for point in torsion_points {
        let frobenius_image = curve.relative_frobenius(&point)?;
        points.push(FrobeniusOnExactTorsionPoint {
            point,
            frobenius_image,
            minimal_absolute_frobenius_fixing_power: None,
        });
    }

    Ok(FrobeniusOnExactTorsionReport {
        exact_order,
        points,
    })
}

/// Computes the absolute Frobenius action `π_p^k` on the exact-`n` rational torsion
/// of a short-Weierstrass curve represented over a finite extension field.
///
/// If the curve coefficients already lie in `F_p` but the points are
/// represented in a larger field such as `F_{p^r}`, then exact-`n` torsion
/// points in `E(F_{p^r})` need not be fixed by the absolute Frobenius `π_p`.
///
/// - enumerate the non-identity rational points of exact order `n`
/// - apply `π_p^k` coordinatewise to each point
/// - record which torsion points are fixed and which are moved
/// - compute, for each point, the smallest positive `d` such that `π_p^d(P)=P`
///
/// Since the current Frobenius API does not yet expose a generic
/// curve-model trait for absolute Frobenius, this function is intentionally
/// specialized to [`ShortWeierstrassCurve<F>`].
///
/// Complexity:  If `N = #E(F_{p^r})`, `m = #E(F_{p^r})[n]_{exact}`, and
/// `s = k mod r`, this funcion performs:
/// - one exhaustive exact-order torsion scan across the rational point set
/// - `m` evaluations of `π_p^k` on points
/// - for each listed point, a scan through powers `1, ..., r` to recover the
///   minimal fixing power inside the represented field
///
/// With the current coordinatewise implementation, each pointwise absolute
/// Frobenius evaluation costs `Θ(s)` repeated field-Frobenius updates, so the
/// direct image computation costs `Θ(ms)`. The additional minimal-fixing-power
/// scan costs `Θ(mr^2)` repeated field-Frobenius updates in the current
/// straightforward implementation, since it evaluates `π_p^d` separately for
/// each `d = 1, ..., r`.
pub fn absolute_frobenius_on_exact_torsion<F: FiniteField + EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    exact_order: usize,
    power: u32,
) -> Result<FrobeniusOnExactTorsionReport<AffinePoint<F>>, CurveError>
where
    F::Elem: Hash,
{
    let torsion_points = points_of_exact_order(curve, exact_order)?;
    let mut points = Vec::with_capacity(torsion_points.len());

    for point in torsion_points {
        let frobenius_image = absolute_frobenius_power_point(curve, &point, power)?;
        let minimal_absolute_frobenius_fixing_power =
            minimal_absolute_frobenius_fixing_power(curve, &point);
        points.push(FrobeniusOnExactTorsionPoint {
            point,
            frobenius_image,
            minimal_absolute_frobenius_fixing_power: Some(minimal_absolute_frobenius_fixing_power),
        });
    }

    Ok(FrobeniusOnExactTorsionReport {
        exact_order,
        points,
    })
}

fn minimal_absolute_frobenius_fixing_power<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
) -> u32 {
    let extension_degree = F::extension_degree().get();
    for power in 1..=extension_degree {
        let image = absolute_frobenius_power_point(curve, point, power)
            .expect("point enumerated from the curve should stay valid under coordinate Frobenius");
        if &image == point {
            return power;
        }
    }

    extension_degree
}
