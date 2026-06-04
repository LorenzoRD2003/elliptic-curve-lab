use num_complex::Complex64;

use super::{ComplexLattice, LatticeIndexPoint};

impl ComplexLattice {
    /// Returns the lattice element `mω₁ + nω₂`.
    pub fn lattice_point(&self, m: i64, n: i64) -> Complex64 {
        self.omega1 * (m as f64) + self.omega2 * (n as f64)
    }

    /// Enumerates lattice points with coefficients in the box
    /// `-radius ≤ m, n ≤ radius`.
    ///
    /// The output order is lexicographic in `(m, n)`, with `m` as the outer
    /// loop and `n` as the inner loop. The origin is included.
    ///
    /// Complexity: `Θ(radius²)` time and `Θ(radius²)` memory, since the full
    /// square box is materialized explicitly.
    pub fn lattice_points_in_box(&self, radius: usize) -> Vec<LatticeIndexPoint> {
        let radius = radius as i64;
        let mut points = Vec::with_capacity(((2 * radius + 1) * (2 * radius + 1)) as usize);

        for m in -radius..=radius {
            for n in -radius..=radius {
                points.push(LatticeIndexPoint {
                    m,
                    n,
                    value: self.lattice_point(m, n),
                });
            }
        }

        points
    }

    /// Enumerates non-zero lattice points with coefficients in the box
    /// `-radius ≤ m, n ≤ radius`.
    ///
    /// This uses the same lexicographic ordering as
    /// [`Self::lattice_points_in_box`] but omits the index pair `(0, 0)`.
    ///
    /// Complexity: `Θ(radius²)` time and `Θ(radius²)` memory, because the
    /// current implementation materializes the full square box first and then
    /// filters out the origin.
    pub fn nonzero_lattice_points_in_box(&self, radius: usize) -> Vec<LatticeIndexPoint> {
        self.lattice_points_in_box(radius)
            .into_iter()
            .filter(|point| point.m != 0 || point.n != 0)
            .collect()
    }
}
