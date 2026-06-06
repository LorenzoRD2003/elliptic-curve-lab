use num_complex::Complex64;

/// Construction error for one Simpson quadrature domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SimpsonQuadratureDomainError {
    NonFiniteInterval,
}

/// One real interval together with the caller-requested composite-Simpson
/// budget.
///
/// This value object stores the raw requested subinterval count and exposes
/// the normalized even budget used by the actual quadrature rule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimpsonQuadratureDomain {
    start: f64,
    end: f64,
    requested_subintervals: usize,
}

impl SimpsonQuadratureDomain {
    /// Builds one validated Simpson quadrature domain.
    ///
    /// Validation currently checks only that both endpoints are finite. The
    /// requested subinterval budget is stored as supplied and normalized later
    /// through [`Self::normalized_subintervals`].
    pub fn new(
        start: f64,
        end: f64,
        requested_subintervals: usize,
    ) -> Result<Self, SimpsonQuadratureDomainError> {
        if !start.is_finite() || !end.is_finite() {
            return Err(SimpsonQuadratureDomainError::NonFiniteInterval);
        }

        Ok(Self {
            start,
            end,
            requested_subintervals,
        })
    }

    /// Returns the left endpoint.
    pub fn start(&self) -> f64 {
        self.start
    }

    /// Returns the right endpoint.
    pub fn end(&self) -> f64 {
        self.end
    }

    /// Returns the raw caller-requested subinterval budget.
    pub fn requested_subintervals(&self) -> usize {
        self.requested_subintervals
    }

    /// Returns the even subinterval budget actually used by composite
    /// Simpson quadrature.
    pub fn normalized_subintervals(&self) -> usize {
        simpson_subinterval_budget(self.requested_subintervals)
    }

    /// Returns the uniform node spacing used after normalization.
    pub fn step_size(&self) -> f64 {
        (self.end - self.start) / self.normalized_subintervals() as f64
    }
}

/// Numerical failure modes for the shared composite-Simpson helper.
#[derive(Clone, Debug, PartialEq)]
pub enum SimpsonIntegrationError<E> {
    Integrand(E),
    NonFiniteIntegrandValue { index: usize, parameter: f64 },
}

/// Normalizes a requested composite-Simpson subinterval count.
///
/// Composite Simpson requires an even number of subintervals. This helper also
/// enforces the smallest meaningful budget, namely two subintervals.
fn simpson_subinterval_budget(requested: usize) -> usize {
    let requested = requested.max(2);
    if requested.is_multiple_of(2) {
        requested
    } else {
        requested + 1
    }
}

/// Integrates one complex-valued function on one validated quadrature domain
/// by composite Simpson quadrature.
///
/// Mathematically, if the normalized subinterval count is `n` and the step is
/// `h = (b-a)/n`, the routine applies the classical formula
/// `∫_a^b f(x) dx ≈ (h / 3) * Σ_{j=0}^n w_j f(a + j h)` with weights
///
/// ```text
/// w_0 = w_n = 1,
/// w_j = 4 for odd interior indices,
/// w_j = 2 for even interior indices.
/// ```
///
/// The closure is evaluated exactly once at each Simpson node, in increasing
/// parameter order, with arguments `(index, parameter)`. That evaluation order
/// is part of the intended contract because some callers, such as
/// branch-continuation routines, need to evolve state sequentially from one
/// sample node to the next.
///
/// Error behavior:
/// - returns [`SimpsonIntegrationError::Integrand`] when the caller-supplied
///   closure fails at one node,
/// - returns [`SimpsonIntegrationError::NonFiniteIntegrandValue`] when the
///   closure succeeds but produces a non-finite complex value.
///
/// Complexity: `Θ(n)`, where `n = domain.normalized_subintervals()`.
pub fn composite_simpson_integrate_complex_in_domain<F, E>(
    domain: &SimpsonQuadratureDomain,
    mut integrand: F,
) -> Result<Complex64, SimpsonIntegrationError<E>>
where
    F: FnMut(usize, f64) -> Result<Complex64, E>,
{
    let subintervals = domain.normalized_subintervals();
    let step = domain.step_size();
    let mut weighted_sum = Complex64::new(0.0, 0.0);

    for index in 0..=subintervals {
        let parameter = domain.start() + step * index as f64;
        let value = integrand(index, parameter).map_err(SimpsonIntegrationError::Integrand)?;
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(SimpsonIntegrationError::NonFiniteIntegrandValue { index, parameter });
        }
        weighted_sum += simpson_weighted_sample(index, subintervals, value);
    }
    Ok(Complex64::new(step / 3.0, 0.0) * weighted_sum)
}

/// Integrates one complex-valued function on one validated quadrature domain
/// by composite Simpson quadrature, using the simpler callback shape
/// `f(parameter)` rather than `f(index, parameter)`.
///
/// Use this wrapper when the caller only depends on the node location and does
/// not need the explicit Simpson node index.
pub fn composite_simpson_integrate_complex_simple_in_domain<F, E>(
    domain: &SimpsonQuadratureDomain,
    mut integrand: F,
) -> Result<Complex64, SimpsonIntegrationError<E>>
where
    F: FnMut(f64) -> Result<Complex64, E>,
{
    composite_simpson_integrate_complex_in_domain(domain, |_, parameter| integrand(parameter))
}

fn simpson_node_weight(index: usize, subintervals: usize) -> f64 {
    if index == 0 || index == subintervals {
        1.0
    } else if index.is_multiple_of(2) {
        2.0
    } else {
        4.0
    }
}

fn simpson_weighted_sample(index: usize, subintervals: usize, value: Complex64) -> Complex64 {
    let weight = simpson_node_weight(index, subintervals);
    Complex64::new(weight, 0.0) * value
}

#[cfg(test)]
mod tests {
    use super::{
        SimpsonIntegrationError, SimpsonQuadratureDomain, SimpsonQuadratureDomainError,
        composite_simpson_integrate_complex_in_domain,
        composite_simpson_integrate_complex_simple_in_domain, simpson_node_weight,
        simpson_subinterval_budget, simpson_weighted_sample,
    };
    use num_complex::Complex64;

    #[test]
    fn simpson_domain_preserves_caller_supplied_fields() {
        let domain = SimpsonQuadratureDomain::new(-1.0, 2.0, 5).unwrap();

        assert_eq!(domain.start(), -1.0);
        assert_eq!(domain.end(), 2.0);
        assert_eq!(domain.requested_subintervals(), 5);
        assert_eq!(domain.normalized_subintervals(), 6);
        assert_eq!(domain.step_size(), 0.5);
    }

    #[test]
    fn simpson_domain_rejects_non_finite_endpoints() {
        let error = SimpsonQuadratureDomain::new(f64::NAN, 1.0, 4).unwrap_err();
        assert_eq!(error, SimpsonQuadratureDomainError::NonFiniteInterval);
    }

    #[test]
    fn simpson_budget_rounds_up_to_an_even_positive_count() {
        assert_eq!(simpson_subinterval_budget(0), 2);
        assert_eq!(simpson_subinterval_budget(1), 2);
        assert_eq!(simpson_subinterval_budget(2), 2);
        assert_eq!(simpson_subinterval_budget(5), 6);
    }

    #[test]
    fn simpson_node_weights_match_the_classical_pattern() {
        assert_eq!(simpson_node_weight(0, 6), 1.0);
        assert_eq!(simpson_node_weight(1, 6), 4.0);
        assert_eq!(simpson_node_weight(2, 6), 2.0);
        assert_eq!(simpson_node_weight(3, 6), 4.0);
        assert_eq!(simpson_node_weight(4, 6), 2.0);
        assert_eq!(simpson_node_weight(5, 6), 4.0);
        assert_eq!(simpson_node_weight(6, 6), 1.0);
    }

    #[test]
    fn weighted_sample_applies_the_expected_simpson_coefficient() {
        let sample = simpson_weighted_sample(3, 6, Complex64::new(1.5, -0.5));
        assert_eq!(sample, Complex64::new(6.0, -2.0));
    }

    #[test]
    fn simpson_integrates_a_constant_exactly() {
        let domain = SimpsonQuadratureDomain::new(0.0, 2.0, 7).unwrap();
        let value = composite_simpson_integrate_complex_in_domain(&domain, |_, _| {
            Ok::<_, ()>(Complex64::new(3.0, -1.0))
        })
        .unwrap();

        assert_eq!(value, Complex64::new(6.0, -2.0));
    }

    #[test]
    fn simpson_simple_domain_wrapper_integrates_a_quadratic_exactly() {
        let domain = SimpsonQuadratureDomain::new(0.0, 2.0, 7).unwrap();
        let value = composite_simpson_integrate_complex_simple_in_domain(&domain, |x| {
            Ok::<_, ()>(Complex64::new((x / 2.0) * (x / 2.0), 0.0))
        })
        .unwrap();

        assert!((value.re - 2.0 / 3.0).abs() <= 1.0e-14);
        assert_eq!(value.im, 0.0);
    }

    #[test]
    fn simpson_calls_the_integrand_in_increasing_node_order() {
        let mut seen = Vec::new();
        let domain = SimpsonQuadratureDomain::new(0.0, 1.0, 5).unwrap();

        let _ = composite_simpson_integrate_complex_in_domain(&domain, |index, x| {
            seen.push((index, x));
            Ok::<_, ()>(Complex64::new(1.0, 0.0))
        })
        .unwrap();

        assert_eq!(seen.len(), 7);
        assert_eq!(seen.first().unwrap().0, 0);
        assert_eq!(seen.last().unwrap().0, 6);
        assert!(seen.windows(2).all(|pair| pair[0].1 < pair[1].1));
    }

    #[test]
    fn simpson_rejects_non_finite_intervals() {
        let error = SimpsonQuadratureDomain::new(f64::NAN, 1.0, 4).unwrap_err();
        assert_eq!(error, SimpsonQuadratureDomainError::NonFiniteInterval);
    }

    #[test]
    fn simpson_propagates_integrand_errors() {
        let domain = SimpsonQuadratureDomain::new(0.0, 1.0, 4).unwrap();
        let error = composite_simpson_integrate_complex_in_domain(&domain, |index, _| {
            if index == 2 {
                Err("boom")
            } else {
                Ok(Complex64::new(0.0, 0.0))
            }
        })
        .unwrap_err();

        assert_eq!(error, SimpsonIntegrationError::Integrand("boom"));
    }

    #[test]
    fn simpson_rejects_non_finite_integrand_values() {
        let domain = SimpsonQuadratureDomain::new(0.0, 1.0, 4).unwrap();
        let error = composite_simpson_integrate_complex_in_domain(&domain, |index, _| {
            if index == 1 {
                Ok::<_, ()>(Complex64::new(f64::INFINITY, 0.0))
            } else {
                Ok::<_, ()>(Complex64::new(0.0, 0.0))
            }
        })
        .unwrap_err();

        assert_eq!(
            error,
            SimpsonIntegrationError::NonFiniteIntegrandValue {
                index: 1,
                parameter: 0.25,
            }
        );
    }
}
