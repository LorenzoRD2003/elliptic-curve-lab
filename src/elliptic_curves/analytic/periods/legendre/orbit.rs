use crate::elliptic_curves::analytic::periods::LegendreParameter;
use num_complex::Complex64;

/// One of the six classical transforms produced by the `S₃` action on a
/// Legendre parameter.
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
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreOrbitElement {
    kind: LegendreOrbitElementKind,
    parameter: LegendreParameter,
}

impl LegendreOrbitElement {
    pub fn kind(&self) -> LegendreOrbitElementKind {
        self.kind
    }

    pub fn parameter(&self) -> &LegendreParameter {
        &self.parameter
    }

    pub fn lambda(&self) -> &Complex64 {
        self.parameter.lambda()
    }
}

/// The full six-element `S₃` orbit of a Legendre parameter.
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

    pub fn elements(&self) -> &[LegendreOrbitElement; 6] {
        &self.elements
    }

    pub fn element(&self, kind: LegendreOrbitElementKind) -> &LegendreOrbitElement {
        &self.elements[kind as usize]
    }

    pub fn values(&self) -> [Complex64; 6] {
        self.elements
            .clone()
            .map(|element| *element.parameter.lambda())
    }
}
