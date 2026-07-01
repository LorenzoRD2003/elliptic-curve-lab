use crate::elliptic_curves::frobenius::FrobeniusTrace;
use crate::isogenies::graphs::IsogenyGraphNodeId;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyFrobeniusRelation {
    degree: usize,
    domain: FrobeniusTrace,
    codomain: FrobeniusTrace,
    same_curve_order: bool,
    same_trace: bool,
}

impl IsogenyFrobeniusRelation {
    pub(crate) fn new(
        degree: usize,
        domain: FrobeniusTrace,
        codomain: FrobeniusTrace,
        same_curve_order: bool,
        same_trace: bool,
    ) -> Self {
        Self {
            degree,
            domain,
            codomain,
            same_curve_order,
            same_trace,
        }
    }

    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn domain(&self) -> &FrobeniusTrace {
        &self.domain
    }
    pub fn codomain(&self) -> &FrobeniusTrace {
        &self.codomain
    }
    pub fn same_curve_order(&self) -> bool {
        self.same_curve_order
    }
    pub fn same_trace(&self) -> bool {
        self.same_trace
    }
    pub fn holds(&self) -> bool {
        self.same_curve_order && self.same_trace
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphNodeFrobeniusData {
    node_id: IsogenyGraphNodeId,
    frobenius_trace: FrobeniusTrace,
    same_curve_order_as_reference: bool,
    same_trace_as_reference: bool,
}

impl IsogenyGraphNodeFrobeniusData {
    pub(crate) fn new(
        node_id: IsogenyGraphNodeId,
        frobenius_trace: FrobeniusTrace,
        same_curve_order_as_reference: bool,
        same_trace_as_reference: bool,
    ) -> Self {
        Self {
            node_id,
            frobenius_trace,
            same_curve_order_as_reference,
            same_trace_as_reference,
        }
    }

    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    pub fn same_curve_order_as_reference(&self) -> bool {
        self.same_curve_order_as_reference
    }

    pub fn same_trace_as_reference(&self) -> bool {
        self.same_trace_as_reference
    }

    pub fn holds(&self) -> bool {
        self.same_curve_order_as_reference && self.same_trace_as_reference
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphFrobeniusReport {
    reference_node: IsogenyGraphNodeId,
    reference: FrobeniusTrace,
    checked_nodes: usize,
    checked_edges: usize,
    node_reports: Vec<IsogenyGraphNodeFrobeniusData>,
    all_same_curve_order: bool,
    all_same_trace: bool,
}

impl IsogenyGraphFrobeniusReport {
    pub(crate) fn new(
        reference_node: IsogenyGraphNodeId,
        reference: FrobeniusTrace,
        checked_nodes: usize,
        checked_edges: usize,
        node_reports: Vec<IsogenyGraphNodeFrobeniusData>,
        all_same_curve_order: bool,
        all_same_trace: bool,
    ) -> Self {
        Self {
            reference_node,
            reference,
            checked_nodes,
            checked_edges,
            node_reports,
            all_same_curve_order,
            all_same_trace,
        }
    }

    pub fn reference_node(&self) -> IsogenyGraphNodeId {
        self.reference_node
    }

    pub fn reference(&self) -> &FrobeniusTrace {
        &self.reference
    }

    pub fn checked_nodes(&self) -> usize {
        self.checked_nodes
    }

    pub fn checked_edges(&self) -> usize {
        self.checked_edges
    }

    pub fn node_reports(&self) -> &[IsogenyGraphNodeFrobeniusData] {
        &self.node_reports
    }

    pub fn all_same_curve_order(&self) -> bool {
        self.all_same_curve_order
    }

    pub fn all_same_trace(&self) -> bool {
        self.all_same_trace
    }

    pub fn holds(&self) -> bool {
        self.all_same_curve_order && self.all_same_trace
    }
}
