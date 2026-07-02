use core::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        torsion::{
            FrobeniusTorsionMatrixError, FrobeniusTorsionMatrixReport, ModNMatrix2, NTorsionBasis,
            TorsionCoordinateTable,
        },
    },
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};

#[cfg(test)]
use crate::elliptic_curves::{
    CurveError,
    frobenius::torsion::{FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport},
    traits::{FiniteGroupCurveModel, RelativeFrobeniusCurveModel},
};

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    pub fn frobenius_matrix_on_n_torsion_basis(
        &self,
        frobenius_trace: FrobeniusTrace,
        basis: NTorsionBasis<AffinePoint<F>>,
    ) -> Result<FrobeniusTorsionMatrixReport<AffinePoint<F>>, FrobeniusTorsionMatrixError>
    where
        F: EnumerableFiniteField + SqrtField,
        F::Elem: PartialEq + Hash,
    {
        let frobenius_power =
            self.validate_frobenius_trace_curve_compatibility(&frobenius_trace)?;
        let coordinate_table = TorsionCoordinateTable::new(self, &basis)?;
        let first_column = self.coordinates_of_frobenius_image_in_basis(
            basis.first(),
            frobenius_power,
            &coordinate_table,
        )?;
        let second_column = self.coordinates_of_frobenius_image_in_basis(
            basis.second(),
            frobenius_power,
            &coordinate_table,
        )?;

        let matrix = ModNMatrix2::from_columns(basis.n(), first_column, second_column)?;

        FrobeniusTorsionMatrixReport::from_matrix_and_trace(frobenius_trace, basis, matrix)
    }

    pub(crate) fn validate_frobenius_trace_curve_compatibility(
        &self,
        frobenius_trace: &FrobeniusTrace,
    ) -> Result<u32, FrobeniusTorsionMatrixError>
    where
        F::Elem: PartialEq,
    {
        let characteristic = F::characteristic().to_biguint();
        if frobenius_trace.base_field().characteristic != characteristic {
            return Err(
                FrobeniusTorsionMatrixError::TraceBaseFieldCharacteristicMismatch {
                    trace_characteristic: frobenius_trace.base_field().characteristic.clone(),
                    curve_characteristic: characteristic.clone(),
                },
            );
        }

        let frobenius_power = frobenius_trace.base_field().extension_degree.get();
        if &self.frobenius_twist_power(frobenius_power).map_err(|_| {
            FrobeniusTorsionMatrixError::FrobeniusTraceDoesNotPreserveCurve {
                extension_degree: frobenius_power,
            }
        })? != self
        {
            return Err(
                FrobeniusTorsionMatrixError::FrobeniusTraceDoesNotPreserveCurve {
                    extension_degree: frobenius_power,
                },
            );
        }

        Ok(frobenius_power)
    }

    pub(crate) fn coordinates_of_frobenius_image_in_basis(
        &self,
        point: &AffinePoint<F>,
        frobenius_power: u32,
        coordinate_table: &TorsionCoordinateTable<AffinePoint<F>>,
    ) -> Result<[usize; 2], FrobeniusTorsionMatrixError>
    where
        F::Elem: PartialEq + Hash,
    {
        let image = self
            .absolute_frobenius_power_point(point, frobenius_power)
            .map_err(|_| FrobeniusTorsionMatrixError::PointNotOnCurve)?;

        coordinate_table
            .coordinates_of(&image)
            .ok_or(FrobeniusTorsionMatrixError::FrobeniusImageOutsideBasisSpan)
    }

    #[cfg(test)]
    pub(crate) fn relative_frobenius_on_exact_torsion(
        &self,
        exact_order: usize,
    ) -> Result<FrobeniusOnExactTorsionReport<AffinePoint<F>>, CurveError>
    where
        F: EnumerableFiniteField + SqrtField,
        F::Elem: Clone + PartialEq,
    {
        let torsion_points = self.points_of_exact_order(exact_order)?;
        let mut points = Vec::with_capacity(torsion_points.len());

        for point in torsion_points {
            let frobenius_image = self.relative_frobenius(&point)?;
            points.push(FrobeniusOnExactTorsionPoint::new(
                point,
                frobenius_image,
                None,
            ));
        }

        Ok(FrobeniusOnExactTorsionReport::new(exact_order, points))
    }

    #[cfg(test)]
    pub(crate) fn absolute_frobenius_on_exact_torsion(
        &self,
        exact_order: usize,
        power: u32,
    ) -> Result<FrobeniusOnExactTorsionReport<AffinePoint<F>>, CurveError>
    where
        F: EnumerableFiniteField + SqrtField,
        F::Elem: Hash,
    {
        let torsion_points = self.points_of_exact_order(exact_order)?;
        let mut points = Vec::with_capacity(torsion_points.len());

        for point in torsion_points {
            let frobenius_image = self.absolute_frobenius_power_point(&point, power)?;
            let minimal_absolute_frobenius_fixing_power =
                self.minimal_absolute_frobenius_fixing_power(&point);
            points.push(FrobeniusOnExactTorsionPoint::new(
                point,
                frobenius_image,
                Some(minimal_absolute_frobenius_fixing_power),
            ));
        }

        Ok(FrobeniusOnExactTorsionReport::new(exact_order, points))
    }

    #[cfg(test)]
    fn minimal_absolute_frobenius_fixing_power(&self, point: &AffinePoint<F>) -> u32 {
        let extension_degree = F::extension_degree().get();
        for power in 1..=extension_degree {
            let image = self.absolute_frobenius_power_point(point, power).expect(
                "point enumerated from the curve should stay valid under coordinate Frobenius",
            );
            if &image == point {
                return power;
            }
        }
        extension_degree
    }
}
