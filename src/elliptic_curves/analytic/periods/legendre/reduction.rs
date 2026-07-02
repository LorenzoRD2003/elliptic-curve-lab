use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    periods::WeierstrassCubicRoots,
    periods::{
        LegendreParameter, LegendreReductionReport,
        legendre::orbit::{LegendreOrbitElementKind, LegendreParameterOrbit},
    },
};
use crate::fields::complex_approx::ComplexApprox;
use crate::numerics::ApproxTolerance;
use num_complex::Complex64;

pub(crate) const LEGENDRE_PERMUTATIONS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// One explicit affine reduction of
/// `4(x - e₁)(x - e₂)(x - e₃)` to Legendre form.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreReduction {
    roots: WeierstrassCubicRoots,
    parameter: LegendreParameter,
    selected_permutation: [usize; 3],
}

impl LegendreReduction {
    /// Builds one deterministic Legendre reduction from an unordered root triple.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        roots.validate_distinct(tolerance)?;

        let selected = roots.choose_legendre_candidate(tolerance)?;
        let parameter = selected.into_parameter(tolerance)?;

        Ok(Self {
            roots: roots.clone(),
            parameter,
            selected_permutation: selected.permutation(),
        })
    }

    pub fn report(&self, tolerance: ApproxTolerance) -> LegendreReductionReport {
        LegendreReductionReport::new(self.clone(), tolerance)
    }

    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    pub fn parameter(&self) -> &LegendreParameter {
        &self.parameter
    }

    pub fn orbit(&self) -> LegendreParameterOrbit {
        self.parameter.orbit()
    }

    pub fn selected_permutation(&self) -> [usize; 3] {
        self.selected_permutation
    }

    pub fn selected_root_triple(&self) -> [&Complex64; 3] {
        let roots = self.roots.roots();
        [
            roots[self.selected_permutation[0]],
            roots[self.selected_permutation[1]],
            roots[self.selected_permutation[2]],
        ]
    }

    pub fn x_translation(&self) -> Complex64 {
        *self.selected_root_triple()[1]
    }

    pub fn x_scale(&self) -> Complex64 {
        *self.selected_root_triple()[0] - *self.selected_root_triple()[1]
    }

    pub fn legendre_rhs_scale_factor(&self) -> Complex64 {
        Complex64::new(4.0, 0.0) * self.x_scale().powu(3)
    }

    pub fn principal_sqrt_x_scale(&self) -> Complex64 {
        self.x_scale().sqrt()
    }

    pub fn legendre_y_scale(&self) -> Complex64 {
        Complex64::new(2.0, 0.0) * self.principal_sqrt_x_scale().powu(3)
    }

    pub fn invariant_differential_scale(&self) -> Complex64 {
        Complex64::new(1.0, 0.0) / (Complex64::new(2.0, 0.0) * self.principal_sqrt_x_scale())
    }

    pub fn legendre_x_from_original_x(&self, original_x: Complex64) -> Complex64 {
        (original_x - self.x_translation()) / self.x_scale()
    }

    pub fn original_x_from_legendre_x(&self, legendre_x: Complex64) -> Complex64 {
        self.x_translation() + self.x_scale() * legendre_x
    }

    pub fn evaluate_legendre_cubic(&self, legendre_x: Complex64) -> Complex64 {
        legendre_x
            * (legendre_x - Complex64::new(1.0, 0.0))
            * (legendre_x - *self.parameter.lambda())
    }

    pub fn evaluate_original_cubic_from_legendre_x(&self, legendre_x: Complex64) -> Complex64 {
        self.legendre_rhs_scale_factor() * self.evaluate_legendre_cubic(legendre_x)
    }

    pub(crate) fn selected_orbit_element_relative_to_input_order(
        &self,
        tolerance: ApproxTolerance,
    ) -> LegendreOrbitElementKind {
        let roots = self.roots().roots();
        let lambda_from_input_order = (roots[2] - roots[1]) / (roots[0] - roots[1]);
        let input_parameter = LegendreParameter::new(lambda_from_input_order)
            .expect("distinct cubic roots must define a finite nonsingular Legendre parameter");

        for element in input_parameter.orbit().elements() {
            if ComplexApprox::eq_with_tolerance(
                element.lambda(),
                self.parameter().lambda(),
                tolerance,
            ) {
                return element.kind();
            }
        }

        unreachable!(
            "selected Legendre parameter must lie in the orbit of the input-order candidate"
        );
    }
}
