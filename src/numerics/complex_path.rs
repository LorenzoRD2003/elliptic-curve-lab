use num_complex::Complex64;

/// One straight line segment in the complex plane.
///
/// This is a small geometric value object for numerical contour routines.
/// The constructor only packages the two endpoints; it does not validate
/// finiteness or non-degeneracy.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexLineSegment {
    start: Complex64,
    end: Complex64,
}

impl ComplexLineSegment {
    /// Builds a line segment from explicit endpoints.
    pub fn new(start: Complex64, end: Complex64) -> Self {
        Self { start, end }
    }

    /// Returns the starting point.
    pub fn start(&self) -> &Complex64 {
        &self.start
    }

    /// Returns the ending point.
    pub fn end(&self) -> &Complex64 {
        &self.end
    }

    /// Returns the displacement vector `end - start`.
    pub fn displacement(&self) -> Complex64 {
        self.end - self.start
    }

    /// Returns the Euclidean length of the segment.
    pub fn length(&self) -> f64 {
        self.displacement().norm()
    }

    /// Returns whether both endpoints are finite complex numbers.
    pub fn is_finite(&self) -> bool {
        self.start.is_finite() && self.end.is_finite()
    }

    /// Returns the affine interpolation point at parameter `t`.
    ///
    /// This uses the standard formula
    /// `γ(t) = (1 - t) start + t end`.
    ///
    /// The method does not restrict `t` to `[0, 1]`; callers may use values
    /// outside that range to continue along the same affine line.
    pub fn point_at(&self, t: f64) -> Complex64 {
        self.start + self.displacement() * t
    }

    /// Returns `subintervals + 1` uniformly spaced sample points including
    /// both endpoints.
    ///
    /// When `subintervals = 0`, this returns the singleton list `[start]`.
    pub fn sample_uniform(&self, subintervals: usize) -> Vec<Complex64> {
        if subintervals == 0 {
            return vec![self.start];
        }

        (0..=subintervals)
            .map(|index| self.point_at(index as f64 / subintervals as f64))
            .collect()
    }
}

/// One ray in the complex plane, described by an origin and a direction angle.
///
/// The direction is `exp(i angle)`, so the ray is
/// `origin + r exp(i angle)` for `r >= 0`.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexRay {
    origin: Complex64,
    angle_radians: f64,
}

impl ComplexRay {
    /// Builds a ray from an origin and an explicit direction angle.
    pub fn new(origin: Complex64, angle_radians: f64) -> Self {
        Self {
            origin,
            angle_radians,
        }
    }

    /// Returns the ray origin.
    pub fn origin(&self) -> &Complex64 {
        &self.origin
    }

    /// Returns the stored direction angle in radians.
    pub fn angle_radians(&self) -> f64 {
        self.angle_radians
    }

    /// Returns the unit complex direction `exp(i angle)`.
    pub fn direction(&self) -> Complex64 {
        Complex64::new(self.angle_radians.cos(), self.angle_radians.sin())
    }

    /// Returns whether the ray origin and angle are finite.
    pub fn is_finite(&self) -> bool {
        self.origin.is_finite() && self.angle_radians.is_finite()
    }

    /// Returns the point at nonnegative distance `r` along the ray.
    pub fn point_at_distance(&self, radius: f64) -> Complex64 {
        self.origin + self.direction() * radius
    }

    /// Returns the point under the standard compactifying map
    /// `r = s / (1 - s)` from `[0, 1)` to `[0, +∞)`.
    ///
    /// In other words, this evaluates
    /// `γ(s) = origin + exp(i angle) * s / (1 - s)`.
    ///
    /// The method does not validate that `s` lies in `[0, 1)`. Values outside
    /// that range simply apply the same rational formula.
    pub fn point_at_compact_parameter(&self, s: f64) -> Complex64 {
        self.point_at_distance(s / (1.0 - s))
    }

    /// Returns the derivative of the compactified parameterization
    /// `γ(s) = origin + exp(i angle) * s / (1 - s)`.
    ///
    /// The derivative is
    /// `γ'(s) = exp(i angle) / (1 - s)^2`.
    ///
    /// As with [`Self::point_at_compact_parameter`], this method does not
    /// validate that `s` lies in `[0, 1)`. It simply evaluates the same
    /// rational formula wherever it is defined.
    pub fn compact_parameter_derivative(&self, s: f64) -> Complex64 {
        self.direction() / (1.0 - s).powi(2)
    }

    /// Returns the compact parameter `s` corresponding to a radial distance
    /// `r` under the standard compactifying map `r = s / (1 - s)`.
    ///
    /// In other words, this evaluates the inverse formula
    /// `s = r / (1 + r)`.
    ///
    /// The method does not validate that `r` is nonnegative. It simply
    /// applies the same rational formula wherever it is defined.
    pub fn compact_parameter_from_distance(&self, radius: f64) -> f64 {
        radius / (1.0 + radius)
    }

    /// Returns `subintervals + 1` sample points obtained by uniformly
    /// subdividing the compact parameter interval `[0, s_max]`.
    ///
    /// This is useful when one wants a finite prefix of the ray while still
    /// using the compactifying parameterization toward infinity.
    ///
    /// When `subintervals = 0`, this returns the singleton list
    /// `[γ(0)] = [origin]`.
    pub fn sample_compact_parameter(&self, s_max: f64, subintervals: usize) -> Vec<Complex64> {
        if subintervals == 0 {
            return vec![self.origin];
        }

        (0..=subintervals)
            .map(|index| {
                let s = s_max * index as f64 / subintervals as f64;
                self.point_at_compact_parameter(s)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::numerics::{ComplexLineSegment, ComplexRay};
    use num_complex::Complex64;

    #[test]
    fn line_segment_preserves_endpoints_and_displacement() {
        let segment = ComplexLineSegment::new(Complex64::new(1.0, -2.0), Complex64::new(4.0, 3.0));

        assert_eq!(segment.start(), &Complex64::new(1.0, -2.0));
        assert_eq!(segment.end(), &Complex64::new(4.0, 3.0));
        assert_eq!(segment.displacement(), Complex64::new(3.0, 5.0));
    }

    #[test]
    fn line_segment_affine_interpolation_hits_midpoint() {
        let segment = ComplexLineSegment::new(Complex64::new(0.0, 0.0), Complex64::new(2.0, 4.0));

        assert_eq!(segment.clone().point_at(0.5), Complex64::new(1.0, 2.0));
        assert_eq!(
            segment.sample_uniform(2),
            vec![
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 2.0),
                Complex64::new(2.0, 4.0),
            ]
        );
    }

    #[test]
    fn line_segment_sampling_with_zero_subintervals_returns_only_the_start() {
        let segment = ComplexLineSegment::new(Complex64::new(2.0, 1.0), Complex64::new(5.0, -1.0));

        assert_eq!(segment.sample_uniform(0), vec![Complex64::new(2.0, 1.0)]);
    }

    #[test]
    fn ray_direction_and_distance_parameterization_match_the_expected_axis() {
        let ray = ComplexRay::new(Complex64::new(1.0, 2.0), 0.0);

        assert_eq!(ray.direction(), Complex64::new(1.0, 0.0));
        assert_eq!(ray.point_at_distance(3.5), Complex64::new(4.5, 2.0));
    }

    #[test]
    fn ray_compact_parameterization_uses_s_over_one_minus_s() {
        let ray = ComplexRay::new(Complex64::new(0.0, 0.0), std::f64::consts::FRAC_PI_2);
        let point = ray.point_at_compact_parameter(0.5);

        assert!((point - Complex64::new(0.0, 1.0)).norm() <= 1.0e-12);
    }

    #[test]
    fn ray_compact_parameter_derivative_matches_the_expected_formula() {
        let ray = ComplexRay::new(Complex64::new(0.0, 0.0), 0.0);
        let derivative = ray.compact_parameter_derivative(0.5);

        assert!((derivative - Complex64::new(4.0, 0.0)).norm() <= 1.0e-12);
    }

    #[test]
    fn ray_compact_parameter_inverse_matches_r_over_one_plus_r() {
        let ray = ComplexRay::new(Complex64::new(0.0, 0.0), 0.0);
        let s = ray.compact_parameter_from_distance(3.0);

        assert!((s - 0.75).abs() <= 1.0e-12);
    }

    #[test]
    fn ray_compact_sampling_with_zero_subintervals_returns_only_the_origin() {
        let ray = ComplexRay::new(Complex64::new(-1.0, 3.0), std::f64::consts::PI);

        assert_eq!(
            ray.sample_compact_parameter(0.75, 0),
            vec![Complex64::new(-1.0, 3.0)]
        );
    }

    #[test]
    fn finite_checks_detect_non_finite_inputs() {
        let finite_segment =
            ComplexLineSegment::new(Complex64::new(0.0, 0.0), Complex64::new(1.0, 1.0));
        let non_finite_segment =
            ComplexLineSegment::new(Complex64::new(f64::NAN, 0.0), Complex64::new(1.0, 1.0));
        let finite_ray = ComplexRay::new(Complex64::new(0.0, 0.0), 1.0);
        let non_finite_ray = ComplexRay::new(Complex64::new(0.0, f64::INFINITY), 1.0);

        assert!(finite_segment.is_finite());
        assert!(!non_finite_segment.is_finite());
        assert!(finite_ray.is_finite());
        assert!(!non_finite_ray.is_finite());
    }
}
