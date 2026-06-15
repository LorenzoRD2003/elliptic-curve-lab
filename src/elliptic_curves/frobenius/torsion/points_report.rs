use std::collections::HashMap;
use std::hash::Hash;

use crate::elliptic_curves::frobenius::orbit::{
    FrobeniusOrbit, orbit_from_successor_by_key, partition_point_orbits_by_key,
};

#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct FrobeniusOnExactTorsionPoint<P> {
    point: P,
    frobenius_image: P,
    minimal_absolute_frobenius_fixing_power: Option<u32>,
}

#[cfg_attr(not(test), allow(dead_code))]
impl<P: PartialEq> FrobeniusOnExactTorsionPoint<P> {
    #[cfg(test)]
    pub fn new(
        point: P,
        frobenius_image: P,
        minimal_absolute_frobenius_fixing_power: Option<u32>,
    ) -> Self {
        Self {
            point,
            frobenius_image,
            minimal_absolute_frobenius_fixing_power,
        }
    }

    pub fn point(&self) -> &P {
        &self.point
    }

    pub fn frobenius_image(&self) -> &P {
        &self.frobenius_image
    }

    pub fn fixed_by_frobenius(&self) -> bool {
        self.point == self.frobenius_image
    }

    pub fn minimal_absolute_frobenius_fixing_power(&self) -> Option<u32> {
        self.minimal_absolute_frobenius_fixing_power
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct FrobeniusOnExactTorsionReport<P> {
    exact_order: usize,
    points: Vec<FrobeniusOnExactTorsionPoint<P>>,
}

#[cfg_attr(not(test), allow(dead_code))]
impl<P: PartialEq> FrobeniusOnExactTorsionReport<P> {
    #[cfg(test)]
    pub fn new(exact_order: usize, points: Vec<FrobeniusOnExactTorsionPoint<P>>) -> Self {
        Self {
            exact_order,
            points,
        }
    }

    pub fn exact_order(&self) -> usize {
        self.exact_order
    }

    pub fn points(&self) -> &[FrobeniusOnExactTorsionPoint<P>] {
        &self.points
    }

    pub fn all_fixed(&self) -> bool {
        self.points.iter().all(|point| point.fixed_by_frobenius())
    }

    pub fn fixed_count(&self) -> usize {
        self.points
            .iter()
            .filter(|point| point.fixed_by_frobenius())
            .count()
    }

    pub fn moved_count(&self) -> usize {
        self.points.len() - self.fixed_count()
    }

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
    pub fn orbit_count(&self) -> usize
    where
        P: Clone + Eq + Hash,
    {
        self.orbits().len()
    }
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
