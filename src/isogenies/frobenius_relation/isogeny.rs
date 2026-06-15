use super::IsogenyFrobeniusRelation;
use crate::elliptic_curves::{CurveError, traits::FrobeniusTraceCurveModel};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::traits::Isogeny;

/// Capability trait for isogenies whose domain and codomain can be compared
/// through their Frobenius-trace packages.
pub trait FrobeniusComparableIsogeny<Domain, Codomain>: Isogeny<Domain, Codomain>
where
    Domain: FrobeniusTraceCurveModel,
    Domain::BaseField:
        EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem> + FiniteField,
    Domain::Point: PartialEq,
    Codomain: FrobeniusTraceCurveModel,
    Codomain::BaseField: EnumerableFiniteField<Elem = Codomain::Elem>
        + SqrtField<Elem = Codomain::Elem>
        + FiniteField,
    Codomain::Point: PartialEq,
{
    /// Compares the domain and codomain Frobenius packages of this isogeny.
    fn frobenius_relation_report(&self) -> Result<IsogenyFrobeniusRelation, CurveError> {
        let domain = self.domain().frobenius_trace()?;
        let codomain = self.codomain().frobenius_trace()?;

        if domain.base_field() != codomain.base_field() {
            return Err(CurveError::IncompatibleFrobeniusIsogenyBaseFields {
                domain_characteristic: domain.base_field().characteristic,
                domain_extension_degree: domain.base_field().extension_degree.get(),
                codomain_characteristic: codomain.base_field().characteristic,
                codomain_extension_degree: codomain.base_field().extension_degree.get(),
            });
        }

        let same_curve_order = domain.curve_order() == codomain.curve_order();
        let same_trace = domain.trace() == codomain.trace();
        Ok(IsogenyFrobeniusRelation::new(
            self.degree(),
            domain,
            codomain,
            same_curve_order,
            same_trace,
        ))
    }
}

impl<I, Domain, Codomain> FrobeniusComparableIsogeny<Domain, Codomain> for I
where
    I: Isogeny<Domain, Codomain>,
    Domain: FrobeniusTraceCurveModel,
    Domain::BaseField:
        EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem> + FiniteField,
    Domain::Point: PartialEq,
    Codomain: FrobeniusTraceCurveModel,
    Codomain::BaseField: EnumerableFiniteField<Elem = Codomain::Elem>
        + SqrtField<Elem = Codomain::Elem>
        + FiniteField,
    Codomain::Point: PartialEq,
{
}
