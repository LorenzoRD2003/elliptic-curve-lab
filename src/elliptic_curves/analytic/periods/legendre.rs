use num_complex::Complex64;

use super::super::AnalyticCurveError;
use super::WeierstrassCubicRoots;
use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

const LEGENDRE_PERMUTATIONS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// One of the six classical transforms produced by the `S₃` action on a
/// Legendre parameter.
///
/// If one ordered root triple `(e₁, e₂, e₃)` defines
///
/// `λ = (e₃ - e₂) / (e₁ - e₂)`,
///
/// then permuting the three roots changes that value by one of six Möbius
/// transforms. For example:
///
/// - swapping `e₁` and `e₃` sends `λ` to `1 / λ`
/// - swapping `e₁` and `e₂` sends `λ` to `1 - λ`
/// - swapping `e₂` and `e₃` sends `λ` to `λ / (λ - 1)`
///
/// The remaining three transforms come from composing those basic
/// transpositions. Generic `λ` values therefore have a six-element orbit,
/// though special symmetric values can collapse several labels to the same
/// complex number.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendreOrbitElementKind {
    /// `λ`
    Lambda,
    /// `1 - λ`
    OneMinusLambda,
    /// `1 / λ`
    ReciprocalLambda,
    /// `1 / (1 - λ)`
    ReciprocalOneMinusLambda,
    /// `(λ - 1) / λ`
    LambdaMinusOneOverLambda,
    /// `λ / (λ - 1)`
    LambdaOverLambdaMinusOne,
}

/// One labeled member of the six-element Legendre orbit.
///
/// The label records which formal `S₃` transform produced the value, not just
/// the resulting complex number. This matters because special parameters may
/// satisfy identities such as `λ = 1 - λ`, so two different orbit labels can
/// coincide numerically.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreOrbitElement {
    kind: LegendreOrbitElementKind,
    parameter: LegendreParameter,
}

impl LegendreOrbitElement {
    /// Returns which classical transform produced this orbit element.
    pub fn kind(&self) -> LegendreOrbitElementKind {
        self.kind
    }

    /// Returns the corresponding Legendre parameter.
    pub fn parameter(&self) -> &LegendreParameter {
        &self.parameter
    }

    /// Returns the corresponding complex number `λ`.
    pub fn lambda(&self) -> &Complex64 {
        self.parameter.lambda()
    }
}

/// The full six-element `S₃` orbit of a Legendre parameter.
///
/// The orbit is stored in the classical fixed order
/// `{λ, 1-λ, 1/λ, 1/(1-λ), (λ-1)/λ, λ/(λ-1)}`.
///
/// Conceptually, this is the orbit of the cross-ratio-like coordinate
/// attached to one ordered triple of roots. Two ordered triples that differ by
/// a permutation need not give the same `λ`, but they always land in this same
/// six-label family.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreParameterOrbit {
    elements: [LegendreOrbitElement; 6],
}

impl LegendreParameterOrbit {
    /// Builds the full orbit from one nonsingular Legendre parameter.
    ///
    /// Complexity: `Θ(1)`.
    pub fn from_parameter(parameter: &LegendreParameter) -> Self {
        let lambda = *parameter.lambda();
        let one = Complex64::new(1.0, 0.0);
        let transforms = [
            (LegendreOrbitElementKind::Lambda, lambda),
            (LegendreOrbitElementKind::OneMinusLambda, one - lambda),
            (LegendreOrbitElementKind::ReciprocalLambda, one / lambda),
            (
                LegendreOrbitElementKind::ReciprocalOneMinusLambda,
                one / (one - lambda),
            ),
            (
                LegendreOrbitElementKind::LambdaMinusOneOverLambda,
                (lambda - one) / lambda,
            ),
            (
                LegendreOrbitElementKind::LambdaOverLambdaMinusOne,
                lambda / (lambda - one),
            ),
        ];

        let elements = transforms.map(|(kind, value)| LegendreOrbitElement {
            kind,
            parameter: LegendreParameter::new(value)
                .expect("Legendre orbit transforms must preserve nonsingularity"),
        });

        Self { elements }
    }

    /// Returns the six orbit elements in their fixed classical order.
    pub fn elements(&self) -> &[LegendreOrbitElement; 6] {
        &self.elements
    }

    /// Returns the element with the requested label.
    pub fn element(&self, kind: LegendreOrbitElementKind) -> &LegendreOrbitElement {
        &self.elements[kind as usize]
    }

    /// Returns just the six complex `λ` values in classical order.
    pub fn values(&self) -> [Complex64; 6] {
        self.elements
            .clone()
            .map(|element| *element.parameter.lambda())
    }
}

/// One chosen Legendre parameter `λ` for a cubic normalized to
/// `y² = x(x - 1)(x - λ)`.
///
/// For an unordered root triple `{e₁, e₂, e₃}`, the quantity
/// `λ = (e₃ - e₂) / (e₁ - e₂)` is not unique. Choosing a different ordering of
/// the same three roots replaces `λ` by one of the six transforms in the
/// `S₃` orbit `{λ, 1 - λ, 1/λ, 1/(1 - λ), (λ - 1)/λ, λ/(λ - 1)}`.
///
/// So `LegendreParameter` should be read as “one chosen representative of a
/// Legendre class”, not as a canonical invariant attached directly to an
/// unordered cubic-root set.
///
/// This type stores one finite, nonsingular representative. The constructor
/// [`Self::from_roots`] chooses that representative deterministically by
/// scanning the six permutation-induced candidates, preferring the one farthest
/// from the singular Legendre locus `{0, 1, ∞}`.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreParameter {
    lambda: Complex64,
}

impl LegendreParameter {
    /// Builds a finite Legendre parameter from an already chosen `λ`.
    ///
    /// This constructor rejects the exact singular values `0` and `1`, and
    /// also rejects non-finite complex values.
    ///
    /// Complexity: `Θ(1)`.
    pub fn new(lambda: Complex64) -> Result<Self, AnalyticCurveError> {
        if !lambda.is_finite()
            || lambda == Complex64::new(0.0, 0.0)
            || lambda == Complex64::new(1.0, 0.0)
        {
            return Err(AnalyticCurveError::InvalidLegendreModulus);
        }

        Ok(Self { lambda })
    }

    /// Builds one deterministic Legendre parameter from an unordered cubic
    /// root triple.
    ///
    /// If the cubic roots are `e₁, e₂, e₃`, any permutation induces a valid
    /// candidate `λ = (e₃ - e₂) / (e₁ - e₂)`.
    ///
    /// Because [`WeierstrassCubicRoots`] intentionally preserves caller order
    /// without claiming that it is canonical, this constructor reorders the
    /// roots internally. It evaluates the six permutation-induced candidates
    /// and chooses a deterministic representative that stays as far as
    /// possible from the singular Legendre set `{0, 1, ∞}` under the score
    ///
    /// `min(|λ|, |1 - λ|, 1 / |λ|)`.
    ///
    /// The final tolerance gate rejects candidates that remain numerically too
    /// close to `0`, `1`, or `∞`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        roots.validate_distinct(tolerance)?;
        let selected = choose_legendre_candidate_from_roots(roots, tolerance)?;
        build_parameter_from_candidate(selected, tolerance)
    }

    /// Returns the stored Legendre parameter `λ`.
    pub fn lambda(&self) -> &Complex64 {
        &self.lambda
    }

    /// Returns the full six-element Legendre orbit of this parameter.
    pub fn orbit(&self) -> LegendreParameterOrbit {
        LegendreParameterOrbit::from_parameter(self)
    }

    /// Returns `1 - λ`.
    pub fn one_minus_lambda(&self) -> Complex64 {
        Complex64::new(1.0, 0.0) - self.lambda
    }

    /// Returns whether `λ ≈ 0` under `tolerance`.
    pub fn is_near_zero(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_zero()
    }

    /// Returns whether `λ ≈ 1` under `tolerance`.
    pub fn is_near_one(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_one()
    }

    /// Returns whether this Legendre parameter is numerically close to the
    /// singular Legendre locus `{0, 1, ∞}`.
    ///
    /// The `∞` branch is detected projectively through a large `|λ|`, using a
    /// reciprocal threshold derived from `tolerance`.
    pub fn is_near_singular(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_singular()
    }
}

/// Coarse numerical conditioning class for one chosen Legendre parameter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendreParameterConditioning {
    /// The chosen `λ` is not close to `0`, `1`, or `∞` under the tolerance.
    Generic,
    /// The chosen `λ` is numerically close to `0`.
    NearZero,
    /// The chosen `λ` is numerically close to `1`.
    NearOne,
    /// The chosen `λ` is numerically close to `∞`, detected through `1/λ`.
    NearInfinity,
}

impl LegendreParameterConditioning {
    /// Returns whether this conditioning class lies near the singular
    /// Legendre locus `{0, 1, ∞}`.
    pub fn is_near_singular(self) -> bool {
        self != Self::Generic
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct LegendreSingularityDiagnostics {
    conditioning: LegendreParameterConditioning,
    singularity_distance: f64,
    is_near_zero: bool,
    is_near_one: bool,
}

impl LegendreSingularityDiagnostics {
    fn analyze(parameter: &LegendreParameter, tolerance: ApproxTolerance) -> Self {
        let is_near_zero = ComplexApprox::eq_with_tolerance(
            parameter.lambda(),
            &Complex64::new(0.0, 0.0),
            tolerance,
        );
        let is_near_one = ComplexApprox::eq_with_tolerance(
            parameter.lambda(),
            &Complex64::new(1.0, 0.0),
            tolerance,
        );
        let is_near_infinity =
            parameter.lambda().norm() >= reciprocal_singularity_threshold(tolerance);
        let conditioning = if is_near_zero {
            LegendreParameterConditioning::NearZero
        } else if is_near_one {
            LegendreParameterConditioning::NearOne
        } else if is_near_infinity {
            LegendreParameterConditioning::NearInfinity
        } else {
            LegendreParameterConditioning::Generic
        };

        Self {
            conditioning,
            singularity_distance: legendre_singularity_distance(parameter.lambda()),
            is_near_zero,
            is_near_one,
        }
    }

    fn conditioning(self) -> LegendreParameterConditioning {
        self.conditioning
    }

    fn singularity_distance(self) -> f64 {
        self.singularity_distance
    }

    fn is_near_zero(self) -> bool {
        self.is_near_zero
    }

    fn is_near_one(self) -> bool {
        self.is_near_one
    }

    fn is_near_singular(self) -> bool {
        self.conditioning.is_near_singular()
    }
}

/// Structured report for one Legendre reduction.
///
/// The report reuses the already computed [`LegendreReduction`] and adds two
/// diagnostic layers:
///
/// - which orbit element kind, relative to the caller-supplied root order,
///   produced the chosen representative
/// - how close the chosen parameter lies to the singular locus `{0, 1, ∞}`
///
/// The first item is intentionally not intrinsic to the Legendre class: it is
/// a label computed relative to the original root order supplied by the caller.
/// When special symmetric values make that orbit label non-unique, the report
/// chooses the first matching label in the fixed classical orbit order.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreReductionReport {
    reduction: LegendreReduction,
    selected_orbit_element_relative_to_input_order: LegendreOrbitElementKind,
    conditioning: LegendreParameterConditioning,
    tolerance: ApproxTolerance,
    singularity_distance: f64,
}

impl LegendreReductionReport {
    /// Builds a Legendre reduction report from an already computed reduction.
    pub fn new(reduction: LegendreReduction, tolerance: ApproxTolerance) -> Self {
        let selected_orbit_element =
            selected_orbit_element_kind_relative_to_input_order(&reduction, tolerance);
        let diagnostics = LegendreSingularityDiagnostics::analyze(reduction.parameter(), tolerance);

        Self {
            reduction,
            selected_orbit_element_relative_to_input_order: selected_orbit_element,
            conditioning: diagnostics.conditioning(),
            tolerance,
            singularity_distance: diagnostics.singularity_distance(),
        }
    }

    /// Builds a report directly from a root triple.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        Ok(Self::new(
            LegendreReduction::from_roots(roots, tolerance)?,
            tolerance,
        ))
    }

    /// Returns the underlying reduction.
    pub fn reduction(&self) -> &LegendreReduction {
        &self.reduction
    }

    /// Returns the chosen Legendre parameter.
    pub fn parameter(&self) -> &LegendreParameter {
        self.reduction.parameter()
    }

    /// Returns which orbit transform of the caller-order candidate produced
    /// the chosen representative.
    ///
    /// This label is defined relative to the input root order used when
    /// building the report. It is not a canonical invariant of the Legendre
    /// class by itself.
    pub fn selected_orbit_element_relative_to_input_order(&self) -> LegendreOrbitElementKind {
        self.selected_orbit_element_relative_to_input_order
    }

    /// Returns the coarse numerical conditioning class for the chosen `λ`.
    pub fn conditioning(&self) -> LegendreParameterConditioning {
        self.conditioning
    }

    /// Returns the tolerance used for the report classification.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the score `min(|λ|, |1-λ|, 1/|λ|)` measuring distance from the
    /// singular locus `{0, 1, ∞}`.
    pub fn singularity_distance(&self) -> f64 {
        self.singularity_distance
    }

    /// Returns whether the chosen `λ` is near the singular locus under the
    /// stored conditioning classification.
    pub fn is_near_singular(&self) -> bool {
        self.conditioning.is_near_singular()
    }
}

/// Classifies one chosen Legendre parameter by its distance to the singular
/// locus `{0, 1, ∞}`. Complexity: `Θ(1)`.
pub fn classify_legendre_parameter_conditioning(
    parameter: &LegendreParameter,
    tolerance: ApproxTolerance,
) -> LegendreParameterConditioning {
    LegendreSingularityDiagnostics::analyze(parameter, tolerance).conditioning()
}

/// Builds a structured Legendre reduction report directly from roots.
///
/// Complexity: `Θ(1)`.
pub fn legendre_reduction_report(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> Result<LegendreReductionReport, AnalyticCurveError> {
    LegendreReductionReport::from_roots(roots, tolerance)
}

/// One explicit affine reduction of
/// `4(x - e₁)(x - e₂)(x - e₃)` to Legendre form.
///
/// If `X = (x - e₂) / (e₁ - e₂)`, then
///
/// `4(x - e₁)(x - e₂)(x - e₃) = 4(e₁ - e₂)^3 X(X - 1)(X - λ)`.
///
/// This struct stores the chosen permutation of the roots, the resulting
/// Legendre parameter, and the induced affine `x`-change of variables.
///
/// The branch-independent quantity is the `y²` scale factor `4(e₁ - e₂)^3`.
/// When an actual `y` scale or invariant-differential scale
/// is needed, this API fixes the principal square root `α = sqrt(e₁ - e₂)` and uses
///
/// - `x = e₂ + (e₁ - e₂) X`
/// - `y = 2 α^3 Y`
///
/// so that `Y² = X(X - 1)(X - λ)`.
///
/// Replacing `α` by `-α` flips the signs of `Y`, `ω₁`, and `ω₂`
/// simultaneously, but does not change the underlying Legendre curve or the
/// resulting period ratio `τ`.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreReduction {
    roots: WeierstrassCubicRoots,
    parameter: LegendreParameter,
    selected_permutation: [usize; 3],
}

impl LegendreReduction {
    /// Builds one deterministic Legendre reduction from an unordered root triple.
    ///
    /// Complexity: `Θ(1)`.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        roots.validate_distinct(tolerance)?;

        let selected = choose_legendre_candidate_from_roots(roots, tolerance)?;
        let parameter = build_parameter_from_candidate(selected, tolerance)?;

        Ok(Self {
            roots: roots.clone(),
            parameter,
            selected_permutation: selected.permutation,
        })
    }

    /// Returns the original root triple.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    /// Returns the chosen Legendre parameter.
    pub fn parameter(&self) -> &LegendreParameter {
        &self.parameter
    }

    /// Returns the full Legendre orbit of the chosen parameter.
    pub fn orbit(&self) -> LegendreParameterOrbit {
        self.parameter.orbit()
    }

    /// Returns the chosen root permutation `[i, j, k]` used to interpret the
    /// stored roots as `(e₁, e₂, e₃)`.
    pub fn selected_permutation(&self) -> [usize; 3] {
        self.selected_permutation
    }

    /// Returns the root triple `(e₁, e₂, e₃)` after applying the selected
    /// reduction permutation.
    ///
    /// This is the triple used by this reduction step. It is not a canonical
    /// ordering of the original cubic roots.
    pub fn selected_root_triple(&self) -> [&Complex64; 3] {
        let roots = self.roots.roots();
        [
            roots[self.selected_permutation[0]],
            roots[self.selected_permutation[1]],
            roots[self.selected_permutation[2]],
        ]
    }

    /// Returns the affine translation `x = e₂ + (e₁ - e₂) X`.
    pub fn x_translation(&self) -> Complex64 {
        *self.selected_root_triple()[1]
    }

    /// Returns the affine `x`-scale `e₁ - e₂`.
    pub fn x_scale(&self) -> Complex64 {
        *self.selected_root_triple()[0] - *self.selected_root_triple()[1]
    }

    /// Returns the factor multiplying the right-hand side after the affine
    /// `x`-change of variables:
    ///
    /// `4(x - e₁)(x - e₂)(x - e₃) = scale * X(X - 1)(X - λ)`.
    ///
    /// At this stage `y` is still the original coordinate, so this is a
    /// scale factor for the right-hand side, not yet for `y` itself.
    pub fn legendre_rhs_scale_factor(&self) -> Complex64 {
        Complex64::new(4.0, 0.0) * self.x_scale().powu(3)
    }

    /// Returns the principal square root `α = sqrt(e₁ - e₂)` of the chosen
    /// `x`-scale.
    ///
    /// This is the branch used to define the concrete `y` scale and the
    /// invariant-differential scale. Flipping to `-α` would change only a
    /// global sign in those quantities.
    pub fn principal_sqrt_x_scale(&self) -> Complex64 {
        self.x_scale().sqrt()
    }

    /// Returns the concrete principal-branch `y` scale `2 α^3`, where
    /// `α = sqrt(e₁ - e₂)` is the principal square root.
    ///
    /// With this convention, the change of variables
    ///
    /// - `x = e₂ + (e₁ - e₂) X`
    /// - `y = legendre_y_scale * Y`
    ///
    /// sends the original cubic to `Y² = X(X-1)(X-λ)`.
    ///
    /// Equivalently, `legendre_y_scale()^2 = legendre_rhs_scale_factor()`.
    pub fn legendre_y_scale(&self) -> Complex64 {
        Complex64::new(2.0, 0.0) * self.principal_sqrt_x_scale().powu(3)
    }

    /// Returns the scale factor relating invariant differentials:
    ///
    /// `dx / y = invariant_differential_scale * dX / Y`.
    ///
    /// Under the convention `α = sqrt(e₁ - e₂)`, this equals `1 / (2 α)`.
    pub fn invariant_differential_scale(&self) -> Complex64 {
        Complex64::new(1.0, 0.0) / (Complex64::new(2.0, 0.0) * self.principal_sqrt_x_scale())
    }

    /// Maps one original `x`-coordinate to the Legendre `X`-coordinate.
    pub fn legendre_x_from_original_x(&self, original_x: Complex64) -> Complex64 {
        (original_x - self.x_translation()) / self.x_scale()
    }

    /// Maps one Legendre `X`-coordinate back to the original `x`-coordinate.
    pub fn original_x_from_legendre_x(&self, legendre_x: Complex64) -> Complex64 {
        self.x_translation() + self.x_scale() * legendre_x
    }

    /// Evaluates the normalized Legendre cubic `X(X - 1)(X - λ)`.
    pub fn evaluate_legendre_cubic(&self, legendre_x: Complex64) -> Complex64 {
        legendre_x
            * (legendre_x - Complex64::new(1.0, 0.0))
            * (legendre_x - *self.parameter.lambda())
    }

    /// Evaluates the original cubic factor through the Legendre variable `X`.
    pub fn evaluate_original_cubic_from_legendre_x(&self, legendre_x: Complex64) -> Complex64 {
        self.legendre_rhs_scale_factor() * self.evaluate_legendre_cubic(legendre_x)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct LegendreCandidate {
    permutation: [usize; 3],
    lambda: Complex64,
    singular_distance: f64,
}

fn build_parameter_from_candidate(
    candidate: LegendreCandidate,
    tolerance: ApproxTolerance,
) -> Result<LegendreParameter, AnalyticCurveError> {
    let parameter = LegendreParameter::new(candidate.lambda)?;
    if parameter.is_near_singular(tolerance) {
        return Err(AnalyticCurveError::InvalidLegendreModulus);
    }

    Ok(parameter)
}

fn choose_legendre_candidate_from_roots(
    roots: &WeierstrassCubicRoots,
    tolerance: ApproxTolerance,
) -> Result<LegendreCandidate, AnalyticCurveError> {
    let stored_roots = roots.roots();
    let mut best: Option<LegendreCandidate> = None;

    for permutation in LEGENDRE_PERMUTATIONS {
        let e1 = *stored_roots[permutation[0]];
        let e2 = *stored_roots[permutation[1]];
        let e3 = *stored_roots[permutation[2]];
        let denominator = e1 - e2;

        if ComplexApprox::is_zero_with_tolerance(&denominator, tolerance) {
            continue;
        }

        let lambda = (e3 - e2) / denominator;
        if !lambda.is_finite() {
            continue;
        }

        let candidate = LegendreCandidate {
            permutation,
            lambda,
            singular_distance: legendre_singularity_distance(&lambda),
        };

        if legendre_candidate_is_better(candidate, best, tolerance) {
            best = Some(candidate);
        }
    }

    best.ok_or(AnalyticCurveError::InvalidLegendreModulus)
}

fn legendre_candidate_is_better(
    candidate: LegendreCandidate,
    current_best: Option<LegendreCandidate>,
    tolerance: ApproxTolerance,
) -> bool {
    let Some(best) = current_best else {
        return true;
    };

    if !tolerance.real_close(candidate.singular_distance, best.singular_distance) {
        return candidate.singular_distance > best.singular_distance;
    }

    let candidate_norm = candidate.lambda.norm();
    let best_norm = best.lambda.norm();
    if !tolerance.real_close(candidate_norm, best_norm) {
        return candidate_norm < best_norm;
    }

    if !tolerance.real_close(candidate.lambda.re, best.lambda.re) {
        return candidate.lambda.re < best.lambda.re;
    }

    candidate.lambda.im < best.lambda.im
}

fn legendre_singularity_distance(lambda: &Complex64) -> f64 {
    let norm = lambda.norm();
    if norm == 0.0 {
        return 0.0;
    }

    norm.min((Complex64::new(1.0, 0.0) - lambda).norm())
        .min(1.0 / norm)
}

fn selected_orbit_element_kind_relative_to_input_order(
    reduction: &LegendreReduction,
    tolerance: ApproxTolerance,
) -> LegendreOrbitElementKind {
    let roots = reduction.roots().roots();
    let lambda_from_input_order = (roots[2] - roots[1]) / (roots[0] - roots[1]);
    let input_parameter = LegendreParameter::new(lambda_from_input_order)
        .expect("distinct cubic roots must define a finite nonsingular Legendre parameter");

    for element in input_parameter.orbit().elements() {
        if ComplexApprox::eq_with_tolerance(
            element.lambda(),
            reduction.parameter().lambda(),
            tolerance,
        ) {
            return element.kind();
        }
    }

    unreachable!("selected Legendre parameter must lie in the orbit of the input-order candidate");
}

fn infinity_proximity_scale(tolerance: ApproxTolerance) -> f64 {
    tolerance
        .absolute
        .max(tolerance.relative)
        .max(f64::EPSILON.sqrt())
}

fn reciprocal_singularity_threshold(tolerance: ApproxTolerance) -> f64 {
    1.0 / infinity_proximity_scale(tolerance)
}
