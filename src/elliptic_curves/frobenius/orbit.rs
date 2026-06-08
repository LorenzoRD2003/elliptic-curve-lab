use std::collections::HashSet;
use std::hash::Hash;

use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::frobenius::{absolute_frobenius_power_point, frobenius_twist_power};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{
    CurveModel, EnumerableCurveModel, RelativeFrobeniusCurveModel,
};
use crate::fields::{EnumerableFiniteField, FiniteField, SqrtField};

/// Orbit of a point under a chosen Frobenius action.
///
/// The stored `points` list records the distinct points `P, π(P), π^2(P), ...`
/// in cyclic order, stopping just before the orbit returns to `P`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusOrbit<P> {
    start: P,
    points: Vec<P>,
}

impl<P> FrobeniusOrbit<P> {
    pub(crate) fn from_points(start: P, points: Vec<P>) -> Self {
        Self { start, points }
    }

    /// Returns the chosen starting point of the orbit.
    pub fn start(&self) -> &P {
        &self.start
    }

    /// Returns the distinct orbit points in cyclic order.
    pub fn points(&self) -> &[P] {
        &self.points
    }

    /// Returns the orbit period.
    ///
    /// Complexity: `Θ(1)`.
    pub fn period(&self) -> usize {
        self.points.len()
    }
}

/// Computes the orbit of one rational point under the relative Frobenius `π_q`.
///
/// In the current represented finite-field setting, `π_q` is the identity on
/// `E(F_q)`, so every relative-Frobenius orbit is a singleton. The helper is
/// still useful as the relative counterpart of the nontrivial absolute orbit
/// surface.
///
/// Complexity: The cost is `Θ(1)` plus the cost of one call to
/// [`RelativeFrobeniusCurveModel::relative_frobenius`].
pub fn relative_frobenius_orbit<E: RelativeFrobeniusCurveModel>(
    curve: &E,
    point: &E::Point,
) -> Result<FrobeniusOrbit<E::Point>, CurveError>
where
    E::BaseField: FiniteField,
    E::Point: Clone + PartialEq,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    orbit_from_successor(point.clone(), 1, |current| {
        curve.relative_frobenius(current)
    })
}

/// Partitions the rational point set into orbits of the relative Frobenius `π_q`.
///
/// Since `π_q` is the identity on `E(F_q)` in the current represented
/// setting, this returns singleton orbits. Even so, it gives a consistent
/// orbit-level API next to the absolute Frobenius version.
///
/// Complexity: If `N = #E(F_q)`, this performs one full point enumeration and one
/// relative-Frobenius orbit check per previously unseen point. In the current
/// setting, where every orbit has period `1`, the total work is
/// `Θ(N)` relative-Frobenius evaluations plus `Θ(N^2)` point comparisons.
pub fn relative_frobenius_orbits_on_points<E: EnumerableCurveModel + RelativeFrobeniusCurveModel>(
    curve: &E,
) -> Result<Vec<FrobeniusOrbit<E::Point>>, CurveError>
where
    E::BaseField: EnumerableFiniteField<Elem = E::Elem> + SqrtField<Elem = E::Elem> + FiniteField,
    E::Point: Clone + PartialEq,
{
    partition_point_orbits(curve.points(), |point| {
        relative_frobenius_orbit(curve, point)
    })
}

/// Computes the orbit of one rational point under the absolute Frobenius `π_p^k`.
///
/// This helper is currently specialized to short-Weierstrass curves over a
/// represented finite field `F_{p^r}`. It requires that the chosen Frobenius
/// power preserve the current curve model; for example, curves with
/// coefficients already in `F_p` satisfy this for `k = 1`.
///
/// If the curve is preserved, then the orbit period divides the order of
/// `π_p^k` on the represented field, namely `r / gcd(r, k)`. The current
/// implementation uses that exact bound to close the orbit without a caller-
/// supplied step limit.
pub fn absolute_frobenius_orbit<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    power: u32,
) -> Result<FrobeniusOrbit<AffinePoint<F>>, CurveError>
where
    F::Elem: Hash,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }
    ensure_absolute_frobenius_preserves_curve(curve, power)?;

    let bound = absolute_frobenius_period_bound::<F>(power);
    orbit_from_successor_by_key(
        point.clone(),
        bound as usize,
        |orbit_point| orbit_point.clone(),
        |current| absolute_frobenius_power_point(curve, current, power),
    )
}

/// Partitions the rational point set into orbits of the absolute Frobenius `π_p^k`.
///
/// This is the first nontrivial orbit surface in the current codebase: for a
/// curve defined over `F_p` but enumerated over `F_{p^r}`, these orbits can
/// have periods strictly larger than `1`.
///
/// Complexity:
/// If `N = #E(F_{p^r})` and `b = r / gcd(r, k)`, this helper performs one full
/// point enumeration and one absolute-Frobenius orbit computation per unseen
/// orbit representative. The total Frobenius work is `Θ(N (k mod r))` in the
/// worst case. The orbit-partition bookkeeping now uses hashed point keys, so
/// it performs `Θ(N)` key lookups and insertions under the usual constant-time
/// hash-table model.
pub fn absolute_frobenius_orbits_on_points<F: FiniteField + EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    power: u32,
) -> Result<Vec<FrobeniusOrbit<AffinePoint<F>>>, CurveError>
where
    F::Elem: Hash,
{
    ensure_absolute_frobenius_preserves_curve(curve, power)?;

    partition_point_orbits_by_key(
        curve.points(),
        |orbit_point| orbit_point.clone(),
        |point| absolute_frobenius_orbit(curve, point, power),
    )
}

pub(crate) fn orbit_from_successor<P, S>(
    start: P,
    max_period: usize,
    mut successor: S,
) -> Result<FrobeniusOrbit<P>, CurveError>
where
    P: Clone + PartialEq,
    S: FnMut(&P) -> Result<P, CurveError>,
{
    let mut points = vec![start.clone()];
    let mut current = start.clone();

    for _ in 0..max_period {
        let next = successor(&current)?;
        if next == start {
            break;
        }
        if points.contains(&next) {
            break;
        }
        points.push(next.clone());
        current = next;
    }

    Ok(FrobeniusOrbit::from_points(start, points))
}

pub(crate) fn orbit_from_successor_by_key<P, K, KeyOf, S>(
    start: P,
    max_period: usize,
    key_of: KeyOf,
    mut successor: S,
) -> Result<FrobeniusOrbit<P>, CurveError>
where
    P: Clone,
    K: Eq + Hash,
    KeyOf: Fn(&P) -> K,
    S: FnMut(&P) -> Result<P, CurveError>,
{
    let start_key = key_of(&start);
    let mut seen_keys = HashSet::from([start_key]);
    let mut points = vec![start.clone()];
    let mut current = start;

    for _ in 0..max_period {
        let next = successor(&current)?;
        let next_key = key_of(&next);
        if !seen_keys.insert(next_key) {
            break;
        }
        points.push(next.clone());
        current = next;
    }

    Ok(FrobeniusOrbit::from_points(points[0].clone(), points))
}

pub(crate) fn partition_point_orbits<P, I, B>(
    points: I,
    mut build_orbit: B,
) -> Result<Vec<FrobeniusOrbit<P>>, CurveError>
where
    P: Clone + PartialEq,
    I: IntoIterator<Item = P>,
    B: FnMut(&P) -> Result<FrobeniusOrbit<P>, CurveError>,
{
    let mut visited = Vec::new();
    let mut orbits = Vec::new();

    for point in points {
        if visited.contains(&point) {
            continue;
        }

        let orbit = build_orbit(&point)?;
        visited.extend(orbit.points().iter().cloned());
        orbits.push(orbit);
    }

    Ok(orbits)
}

pub(crate) fn partition_point_orbits_by_key<P, K, I, KeyOf, B>(
    points: I,
    key_of: KeyOf,
    mut build_orbit: B,
) -> Result<Vec<FrobeniusOrbit<P>>, CurveError>
where
    P: Clone,
    K: Eq + Hash,
    I: IntoIterator<Item = P>,
    KeyOf: Fn(&P) -> K,
    B: FnMut(&P) -> Result<FrobeniusOrbit<P>, CurveError>,
{
    let mut visited_keys = HashSet::new();
    let mut orbits = Vec::new();

    for point in points {
        let point_key = key_of(&point);
        if visited_keys.contains(&point_key) {
            continue;
        }

        let orbit = build_orbit(&point)?;
        for orbit_point in orbit.points() {
            visited_keys.insert(key_of(orbit_point));
        }
        orbits.push(orbit);
    }

    Ok(orbits)
}

fn ensure_absolute_frobenius_preserves_curve<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    power: u32,
) -> Result<(), CurveError> {
    let twist = frobenius_twist_power(curve, power)?;
    if F::eq(curve.a(), twist.a()) && F::eq(curve.b(), twist.b()) {
        Ok(())
    } else {
        Err(CurveError::AbsoluteFrobeniusDoesNotPreserveCurve { power })
    }
}

fn absolute_frobenius_period_bound<F: FiniteField>(power: u32) -> u32 {
    let extension_degree = F::extension_degree().get();
    let reduced_power = power % extension_degree;
    if reduced_power == 0 {
        return 1;
    }

    extension_degree / gcd_u32(extension_degree, reduced_power)
}

fn gcd_u32(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }

    left
}
