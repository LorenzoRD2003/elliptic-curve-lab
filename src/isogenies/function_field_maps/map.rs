use core::fmt;

use crate::elliptic_curves::{
    ShortWeierstrassCurve, ShortWeierstrassFunction, ShortWeierstrassFunctionField,
};
use crate::fields::{Field, RationalFunction};
use crate::isogenies::IsogenyError;
use crate::polynomials::DensePolynomial;

/// Pullback map `φ* : F(E') -> F(E)` between short-Weierstrass function fields.
///
/// This type records the function-field data attached to a morphism
///
/// `phi : E -> E'`
///
/// between validated short-Weierstrass curves. Instead of starting from point
/// formulas directly, it stores the images of the codomain coordinate functions:
///
/// `φ*(x') = X_φ ∈ F(E)`, and `φ*(y') = Y_φ ∈ F(E)`
///
/// The constructor validates two basic invariants:
///
/// - both stored pullbacks live on the declared domain curve `E`
/// - the codomain relation is respected after substitution:
///   `Y_φ^2 = X_φ^3 + a' X_φ + b'`
///
/// This is enough to model the induced algebra map on the current educational
/// function-field layer. It does **not** yet certify that the data comes from
/// a genuine finite isogeny, nor that the induced map is injective on
/// function fields.
#[derive(Clone)]
pub struct ShortWeierstrassFunctionFieldMap<F: Field> {
    domain_curve: ShortWeierstrassCurve<F>,
    codomain_curve: ShortWeierstrassCurve<F>,
    x_pullback: ShortWeierstrassFunction<F>,
    y_pullback: ShortWeierstrassFunction<F>,
}

impl<F: Field> ShortWeierstrassFunctionFieldMap<F>
where
    F::Elem: PartialEq,
{
    /// Builds a pullback map `φ* : F(E') -> F(E)` from the images of `x'`
    /// and `y'`.
    ///
    /// The stored `x_pullback` and `y_pullback` must be elements of the
    /// domain function field `F(E)`, and together they must satisfy the
    /// codomain equation after substitution.
    pub fn new(
        domain_curve: ShortWeierstrassCurve<F>,
        codomain_curve: ShortWeierstrassCurve<F>,
        x_pullback: ShortWeierstrassFunction<F>,
        y_pullback: ShortWeierstrassFunction<F>,
    ) -> Result<Self, IsogenyError> {
        Self::ensure_pullbacks_live_on_domain(&domain_curve, &x_pullback, &y_pullback)?;
        Self::ensure_codomain_equation_holds(&codomain_curve, &x_pullback, &y_pullback)?;

        Ok(Self {
            domain_curve,
            codomain_curve,
            x_pullback,
            y_pullback,
        })
    }

    /// Returns the domain curve `E` of the map `φ: E -> E'`.
    pub fn domain_curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.domain_curve
    }

    /// Returns the codomain curve `E` of the map `φ: E -> E'`.
    pub fn codomain_curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.codomain_curve
    }

    /// Returns the stored image `φ*(x')`.
    pub fn x_pullback(&self) -> &ShortWeierstrassFunction<F> {
        &self.x_pullback
    }

    /// Returns the stored image `φ*(y')`.
    pub fn y_pullback(&self) -> &ShortWeierstrassFunction<F> {
        &self.y_pullback
    }

    /// Returns the ambient domain function field `F(E)`.
    pub fn domain_function_field(&self) -> ShortWeierstrassFunctionField<F> {
        ShortWeierstrassFunctionField::new(self.domain_curve.clone())
    }

    /// Returns the ambient codomain function field `F(E')`.
    pub fn codomain_function_field(&self) -> ShortWeierstrassFunctionField<F> {
        ShortWeierstrassFunctionField::new(self.codomain_curve.clone())
    }

    /// Pulls back a polynomial in the codomain `x'`-coordinate.
    ///
    /// If `p(T) = c_0 + c_1 T + ... + c_n T^n`, this returns
    /// `p(φ*(x'))` computed in the domain function field `F(E)`.
    pub fn pullback_polynomial(
        &self,
        polynomial: &DensePolynomial<F>,
    ) -> Result<ShortWeierstrassFunction<F>, IsogenyError> {
        ShortWeierstrassFunction::<F>::evaluate_polynomial_in_x(polynomial, &self.x_pullback)
            .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)
    }

    /// Pulls back a rational function in the codomain coordinate `x'`.
    ///
    /// For `r(x') = p(x') / q(x')`, this computes
    ///
    /// `φ*(r) = p(φ*(x')) / q(φ*(x'))`
    ///
    /// inside `F(E)`.
    pub fn pullback_rational_function(
        &self,
        function: &RationalFunction<F>,
    ) -> Result<ShortWeierstrassFunction<F>, IsogenyError> {
        ShortWeierstrassFunction::<F>::substitute_rational_function_in_x(function, &self.x_pullback)
            .map_err(|error| match error {
                crate::elliptic_curves::CurveError::NonInvertibleFunctionFieldElement => {
                    IsogenyError::FunctionFieldMapDenominatorMapsToZero
                }
                _ => IsogenyError::FunctionFieldMapPullbackCurveMismatch,
            })
    }

    /// Pulls back an element `A(x') + y' B(x')` of the codomain function field.
    ///
    /// Writing the codomain element in the basis `1, y'` over `F(x')`, the
    /// pullback is computed by substitution:
    ///
    /// `φ*(A(x') + y' B(x')) = A(φ*(x')) + φ*(y') * B(φ*(x'))`.
    pub fn pullback_function(
        &self,
        function: &ShortWeierstrassFunction<F>,
    ) -> Result<ShortWeierstrassFunction<F>, IsogenyError> {
        if function.curve() != &self.codomain_curve {
            return Err(IsogenyError::FunctionFieldMapSourceCurveMismatch);
        }

        let pulled_a = self.pullback_rational_function(function.a_part())?;
        let pulled_b = self.pullback_rational_function(function.b_part())?;
        let y_term = self
            .y_pullback
            .mul(&pulled_b)
            .map_err(|_| IsogenyError::FunctionFieldMapSourceCurveMismatch)?;

        pulled_a
            .add(&y_term)
            .map_err(|_| IsogenyError::FunctionFieldMapSourceCurveMismatch)
    }

    /// Composes pullbacks contravariantly.
    ///
    /// If `self` represents `φ* : F(E') -> F(E)` and `next` represents
    /// `Ψ* : F(E'') -> F(E')`, then the returned map represents
    ///
    /// `(Ψ o φ)* = φ* o Ψ* : F(E'') -> F(E)`.
    pub fn compose(&self, next: &Self) -> Result<Self, IsogenyError> {
        if self.codomain_curve != next.domain_curve {
            return Err(IsogenyError::CompositionDomainCodomainMismatch);
        }

        let x_pullback = self.pullback_function(next.x_pullback())?;
        let y_pullback = self.pullback_function(next.y_pullback())?;

        Self::new(
            self.domain_curve.clone(),
            next.codomain_curve.clone(),
            x_pullback,
            y_pullback,
        )
    }

    fn ensure_pullbacks_live_on_domain(
        domain_curve: &ShortWeierstrassCurve<F>,
        x_pullback: &ShortWeierstrassFunction<F>,
        y_pullback: &ShortWeierstrassFunction<F>,
    ) -> Result<(), IsogenyError> {
        if x_pullback.curve() != domain_curve || y_pullback.curve() != domain_curve {
            return Err(IsogenyError::FunctionFieldMapPullbackCurveMismatch);
        }

        Ok(())
    }

    fn ensure_codomain_equation_holds(
        codomain_curve: &ShortWeierstrassCurve<F>,
        x_pullback: &ShortWeierstrassFunction<F>,
        y_pullback: &ShortWeierstrassFunction<F>,
    ) -> Result<(), IsogenyError> {
        let lhs = y_pullback
            .mul(y_pullback)
            .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)?;
        let rhs = codomain_rhs_function(codomain_curve, x_pullback)?;

        if lhs != rhs {
            return Err(IsogenyError::FunctionFieldMapCodomainEquationViolation);
        }

        Ok(())
    }
}

fn codomain_rhs_function<F: Field>(
    codomain_curve: &ShortWeierstrassCurve<F>,
    x_pullback: &ShortWeierstrassFunction<F>,
) -> Result<ShortWeierstrassFunction<F>, IsogenyError> {
    let x_squared = x_pullback
        .mul(x_pullback)
        .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)?;
    let x_cubed = x_squared
        .mul(x_pullback)
        .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)?;
    let a_times_x = x_pullback
        .mul(&ShortWeierstrassFunction::from_rational_function(
            x_pullback.curve().clone(),
            RationalFunction::<F>::constant(codomain_curve.a().clone()),
        ))
        .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)?;
    let constant_term = ShortWeierstrassFunction::from_rational_function(
        x_pullback.curve().clone(),
        RationalFunction::<F>::constant(codomain_curve.b().clone()),
    );

    x_cubed
        .add(&a_times_x)
        .and_then(|partial| partial.add(&constant_term))
        .map_err(|_| IsogenyError::FunctionFieldMapPullbackCurveMismatch)
}

impl<F> PartialEq for ShortWeierstrassFunctionFieldMap<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.domain_curve == other.domain_curve
            && self.codomain_curve == other.codomain_curve
            && self.x_pullback == other.x_pullback
            && self.y_pullback == other.y_pullback
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassFunctionFieldMap<F>
where
    F::Elem: fmt::Debug,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ShortWeierstrassFunctionFieldMap")
            .field("domain_curve", &self.domain_curve)
            .field("codomain_curve", &self.codomain_curve)
            .field("x_pullback", &self.x_pullback)
            .field("y_pullback", &self.y_pullback)
            .finish()
    }
}
