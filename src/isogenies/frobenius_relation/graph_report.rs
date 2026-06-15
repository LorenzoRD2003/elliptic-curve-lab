use crate::elliptic_curves::{CurveError, traits::FrobeniusTraceCurveModel};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphError, IsogenyGraphNodeId,
};

use super::{IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData};

/// Capability trait for isogeny graphs whose node representatives can be
/// compared through Frobenius-trace packages.
pub trait FrobeniusComparableIsogenyGraph<C: FrobeniusTraceCurveModel + GraphCurveModel>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + core::hash::Hash + PartialEq,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    /// Compares the Frobenius packages of all node representatives in the graph.
    fn frobenius_relation_report(&self) -> Result<IsogenyGraphFrobeniusReport, IsogenyGraphError>;
}

impl<C: FrobeniusTraceCurveModel + GraphCurveModel> FrobeniusComparableIsogenyGraph<C>
    for IsogenyGraph<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + core::hash::Hash + PartialEq,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    fn frobenius_relation_report(&self) -> Result<IsogenyGraphFrobeniusReport, IsogenyGraphError> {
        let reference_node = self
            .nodes()
            .first()
            .ok_or(IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(0)))?;
        let reference = reference_node.representative().frobenius_trace()?;

        let mut node_reports = Vec::with_capacity(self.node_count());
        let mut all_same_curve_order = true;
        let mut all_same_trace = true;

        for node in self.nodes() {
            let frobenius_trace = node.representative().frobenius_trace()?;
            if frobenius_trace.base_field() != reference.base_field() {
                return Err(CurveError::IncompatibleFrobeniusIsogenyBaseFields {
                    domain_characteristic: reference.base_field().characteristic,
                    domain_extension_degree: reference.base_field().extension_degree.get(),
                    codomain_characteristic: frobenius_trace.base_field().characteristic,
                    codomain_extension_degree: frobenius_trace.base_field().extension_degree.get(),
                }
                .into());
            }
            let same_curve_order_as_reference =
                frobenius_trace.curve_order() == reference.curve_order();
            let same_trace_as_reference = frobenius_trace.trace() == reference.trace();
            all_same_curve_order &= same_curve_order_as_reference;
            all_same_trace &= same_trace_as_reference;
            node_reports.push(IsogenyGraphNodeFrobeniusData::new(
                node.id(),
                frobenius_trace,
                same_curve_order_as_reference,
                same_trace_as_reference,
            ));
        }

        Ok(IsogenyGraphFrobeniusReport::new(
            reference_node.id(),
            reference,
            self.node_count(),
            self.edge_count(),
            node_reports,
            all_same_curve_order,
            all_same_trace,
        ))
    }
}
