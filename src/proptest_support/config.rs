/// Shared sizing knobs for internal property-testing strategies.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProptestSupportConfig {
    /// Field-side configuration.
    pub fields: FieldStrategyConfig,
    /// Polynomial-side configuration.
    pub polynomials: PolynomialStrategyConfig,
    /// Curve-side configuration.
    pub curves: CurveStrategyConfig,
    /// Isogeny-side configuration.
    pub isogenies: IsogenyStrategyConfig,
    /// Analytic-side configuration.
    pub analytic: AnalyticStrategyConfig,
}

/// Knobs shared by field-element strategies.
#[derive(Clone, Copy, Debug)]
pub struct FieldStrategyConfig {
    /// Inclusive absolute bound for signed integer-like samples.
    pub max_abs_i64: i64,
    /// Inclusive bound for real parts in approximate complex samples.
    pub max_real_norm: f64,
    /// Inclusive bound for imaginary parts in approximate complex samples.
    pub max_imaginary_norm: f64,
}

impl Default for FieldStrategyConfig {
    fn default() -> Self {
        Self {
            max_abs_i64: 8,
            max_real_norm: 3.0,
            max_imaginary_norm: 3.0,
        }
    }
}

/// Knobs shared by polynomial strategies.
#[derive(Clone, Copy, Debug)]
pub struct PolynomialStrategyConfig {
    /// Maximum dense coefficient length.
    pub max_len: usize,
    /// Maximum number of sparse or multivariate terms.
    pub max_terms: usize,
    /// Maximum univariate degree used in sparse generators.
    pub max_degree: usize,
    /// Maximum exponent in multivariate monomials.
    pub max_exponent: usize,
    /// Ambient arity for multivariate polynomials.
    pub arity: usize,
    /// Whether non-zero dense polynomials should keep a non-zero leading
    /// coefficient after generation.
    pub require_nonzero_leading_coefficient: bool,
}

impl Default for PolynomialStrategyConfig {
    fn default() -> Self {
        Self {
            max_len: 6,
            max_terms: 6,
            max_degree: 6,
            max_exponent: 4,
            arity: 3,
            require_nonzero_leading_coefficient: false,
        }
    }
}

/// Knobs shared by elliptic-curve strategies.
#[derive(Clone, Copy, Debug)]
pub struct CurveStrategyConfig {
    /// Whether sampled point cases may use the identity point.
    pub include_identity_points: bool,
    /// Maximum division-polynomial index requested by generic fixtures.
    pub max_division_index: usize,
}

impl Default for CurveStrategyConfig {
    fn default() -> Self {
        Self {
            include_identity_points: true,
            max_division_index: 6,
        }
    }
}

/// Knobs shared by isogeny strategies.
#[derive(Clone, Copy, Debug)]
pub struct IsogenyStrategyConfig {
    /// Preferred non-trivial short-Weierstrass scaling factors over `F41`.
    pub preferred_bridge_scales: [i64; 3],
}

impl Default for IsogenyStrategyConfig {
    fn default() -> Self {
        Self {
            preferred_bridge_scales: [2, 3, 5],
        }
    }
}

/// Knobs shared by analytic strategies.
#[derive(Clone, Copy, Debug)]
pub struct AnalyticStrategyConfig {
    /// Inclusive bound for the real part of `τ`.
    pub max_real_part: f64,
    /// Strict lower bound for the imaginary part of `τ`.
    pub min_imaginary_part: f64,
    /// Inclusive upper bound for the imaginary part of `τ`.
    pub max_imaginary_part: f64,
}

impl Default for AnalyticStrategyConfig {
    fn default() -> Self {
        Self {
            max_real_part: 2.0,
            min_imaginary_part: 0.1,
            max_imaginary_part: 3.0,
        }
    }
}

pub(crate) fn touch_config_inventory() {
    let support = ProptestSupportConfig::default();
    let _ = support.fields.max_abs_i64;
    let _ = support.fields.max_real_norm;
    let _ = support.fields.max_imaginary_norm;
    let _ = support.polynomials.max_len;
    let _ = support.polynomials.max_terms;
    let _ = support.polynomials.max_degree;
    let _ = support.polynomials.max_exponent;
    let _ = support.polynomials.arity;
    let _ = support.polynomials.require_nonzero_leading_coefficient;
    let _ = support.curves.include_identity_points;
    let _ = support.curves.max_division_index;
    let _ = support.isogenies.preferred_bridge_scales;
    let _ = support.analytic.max_real_part;
    let _ = support.analytic.min_imaginary_part;
    let _ = support.analytic.max_imaginary_part;
}
