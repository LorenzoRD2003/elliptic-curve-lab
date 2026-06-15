use std::collections::HashSet;
use std::hash::Hash;

use crate::elliptic_curves::CurveError;

#[cfg(test)]
use crate::elliptic_curves::traits::RelativeFrobeniusCurveModel;
#[cfg(test)]
use crate::fields::traits::FiniteField;

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

#[cfg(test)]
pub(crate) fn relative_frobenius_orbit<E: RelativeFrobeniusCurveModel>(
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

    Ok(FrobeniusOrbit::from_points(
        point.clone(),
        vec![curve.relative_frobenius(point)?],
    ))
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
