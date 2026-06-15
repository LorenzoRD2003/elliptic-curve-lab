use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, FundamentalParallelogramCoordinate,
    torsion::{TorusTorsionIndex, TorusTorsionPoint},
};

impl ComplexLattice {
    /// Enumerates all torus `n`-torsion points in lexicographic `(a, b)` order.
    ///
    /// Relative to `Λ = ℤω₁ + ℤω₂`, these are the `n²` classes represented by
    /// `z = (a/n)ω₁ + (b/n)ω₂ mod Λ`, with `0 ≤ a, b < n`.
    ///
    /// The returned order is lexicographic in `(a, b)`, with `a` as the outer loop
    /// and `b` as the inner loop. This keeps examples and tests stable.
    ///
    /// Complexity: `Θ(n²)` time and `Θ(n²)` memory, since the full `n × n` grid of
    /// reduced classes is materialized explicitly.
    pub fn torus_n_torsion_points(
        &self,
        n: usize,
    ) -> Result<Vec<TorusTorsionPoint>, AnalyticCurveError> {
        if n == 0 {
            return Err(AnalyticCurveError::InvalidTorusTorsionIndex);
        }

        let mut points = Vec::with_capacity(n.saturating_mul(n));

        for a in 0..n {
            for b in 0..n {
                points.push(torus_torsion_point_from_index(
                    self,
                    TorusTorsionIndex::new(a, b, n)?,
                ));
            }
        }

        Ok(points)
    }

    /// Enumerates the primitive torus `n`-torsion points in lexicographic order.
    ///
    /// Here “primitive” means exact torus order `n`, equivalently
    /// `gcd(a, b, n) = 1` for the reduced class `(a, b; n)`.
    ///
    /// Complexity: `Θ(n²)` time and `Θ(n²)` memory, because the full torus
    /// `n`-torsion grid is enumerated first and then filtered.
    pub fn primitive_torus_n_torsion_points(
        &self,
        n: usize,
    ) -> Result<Vec<TorusTorsionPoint>, AnalyticCurveError> {
        Ok(self
            .torus_n_torsion_points(n)?
            .into_iter()
            .filter(|point| point.is_primitive())
            .collect())
    }
}

pub(super) fn torus_torsion_point_from_index(
    lattice: &ComplexLattice,
    index: TorusTorsionIndex,
) -> TorusTorsionPoint {
    let denominator = index.n() as f64;
    let coordinate = FundamentalParallelogramCoordinate::new(
        index.a() as f64 / denominator,
        index.b() as f64 / denominator,
    )
    .expect("validated reduced torus-torsion indices lie in the unit square");
    let z = lattice.point_from_fundamental_coordinates(coordinate.clone());

    TorusTorsionPoint {
        index,
        coordinate,
        z,
    }
}
