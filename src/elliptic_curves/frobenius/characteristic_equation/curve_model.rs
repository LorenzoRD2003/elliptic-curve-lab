use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        FrobeniusCharacteristicPolynomial, FrobeniusTrace,
        characteristic_equation::{
            FrobeniusCharacteristicEquationCheck, FrobeniusCharacteristicEquationExhaustiveReport,
            FrobeniusCharacteristicEquationTerms,
        },
    },
    traits::{FrobeniusTraceCurveModel, RelativeFrobeniusCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};

pub trait FrobeniusCharacteristicEquationCurveModel: RelativeFrobeniusCurveModel
where
    Self::BaseField: FiniteField,
    Self::Point: Clone + PartialEq,
{
    fn compute_characteristic_equation_terms(
        &self,
        point: &Self::Point,
        characteristic_polynomial: &FrobeniusCharacteristicPolynomial,
    ) -> Result<FrobeniusCharacteristicEquationTerms<Self::Point>, CurveError> {
        FrobeniusTrace::assert_base_field_compatible_with_curve::<Self>(
            characteristic_polynomial.base_field(),
        )?;

        let pi_q = self.relative_frobenius(point)?;
        let pi_q_squared = self.relative_frobenius_squared(point)?;
        let trace_term = self.mul_scalar_signed(&pi_q, characteristic_polynomial.trace())?;
        let q_scalar = u64::try_from(characteristic_polynomial.field_order()).map_err(|_| {
            CurveError::UnsupportedFrobeniusFieldOrder {
                field_order: characteristic_polynomial.field_order(),
            }
        })?;
        let q_times_point = self.mul_scalar(point, q_scalar)?;
        let lhs_without_q = self.sub(&pi_q_squared, &trace_term)?;
        let lhs = self.add(&lhs_without_q, &q_times_point)?;

        Ok(FrobeniusCharacteristicEquationTerms::new(
            pi_q,
            pi_q_squared,
            trace_term,
            q_times_point,
            lhs,
        ))
    }

    fn verify_frobenius_characteristic_equation_at_point(
        &self,
        point: &Self::Point,
        characteristic_polynomial: &FrobeniusCharacteristicPolynomial,
    ) -> Result<FrobeniusCharacteristicEquationCheck<Self::Point>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let terms = self.compute_characteristic_equation_terms(point, characteristic_polynomial)?;
        let holds = self.is_identity(terms.lhs());
        Ok(FrobeniusCharacteristicEquationCheck::from_terms(
            point.clone(),
            terms,
            holds,
        ))
    }

    fn verify_frobenius_characteristic_equation_exhaustive(
        &self,
    ) -> Result<FrobeniusCharacteristicEquationExhaustiveReport<Self::Point>, CurveError>
    where
        Self: FrobeniusTraceCurveModel,
        Self::BaseField:
            EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem> + FiniteField,
    {
        let frobenius_trace = self.frobenius_trace()?;
        frobenius_trace.assert_compatible_with_curve::<Self>()?;
        let points = self.points();
        let checked_points = points.len();
        let characteristic_polynomial = frobenius_trace.characteristic_polynomial();
        let mut failed_checks = Vec::new();

        for point in points {
            let check = self.verify_frobenius_characteristic_equation_at_point(
                &point,
                &characteristic_polynomial,
            )?;
            if !check.holds() {
                failed_checks.push(check);
            }
        }

        Ok(FrobeniusCharacteristicEquationExhaustiveReport::new(
            frobenius_trace,
            checked_points,
            failed_checks,
        ))
    }
}

impl<E: RelativeFrobeniusCurveModel> FrobeniusCharacteristicEquationCurveModel for E
where
    E::BaseField: FiniteField,
    E::Point: Clone + PartialEq,
{
}
