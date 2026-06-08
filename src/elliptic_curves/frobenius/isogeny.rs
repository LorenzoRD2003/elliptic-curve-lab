use crate::elliptic_curves::frobenius::FrobeniusTrace;
use crate::elliptic_curves::{CurveError, FrobeniusTraceCurveModel};
use crate::fields::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::Isogeny;
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphError, IsogenyGraphNodeId,
};

/// Frobenius relation between the domain and codomain of an isogeny over `F_q`.
///
/// If two elliptic curves over the same finite field `F_q` are isogenous over
/// `F_q`, then they have the same number of `F_q`-rational points. Equivalently,
/// their Frobenius traces coincide.
///
/// In the current educational implementation, this report is derived by
/// exhaustively counting rational points on the domain and codomain curves and
/// comparing the resulting Frobenius traces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyFrobeniusRelation {
    degree: usize,
    domain: FrobeniusTrace,
    codomain: FrobeniusTrace,
    same_curve_order: bool,
    same_trace: bool,
}

impl IsogenyFrobeniusRelation {
    /// Returns the degree of the isogeny.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the Frobenius trace package of the domain curve.
    pub fn domain(&self) -> &FrobeniusTrace {
        &self.domain
    }

    /// Returns the Frobenius trace package of the codomain curve.
    pub fn codomain(&self) -> &FrobeniusTrace {
        &self.codomain
    }

    /// Returns whether the two curves have the same order over `F_q`.
    pub fn same_curve_order(&self) -> bool {
        self.same_curve_order
    }

    /// Returns whether the two Frobenius traces coincide.
    pub fn same_trace(&self) -> bool {
        self.same_trace
    }

    /// Returns whether both the order and trace comparisons hold.
    ///
    /// Complexity: `Θ(1)`.
    pub fn holds(&self) -> bool {
        self.same_curve_order && self.same_trace
    }
}

/// Frobenius data attached to one representative node of an isogeny graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphNodeFrobeniusData {
    node_id: IsogenyGraphNodeId,
    frobenius_trace: FrobeniusTrace,
    same_curve_order_as_reference: bool,
    same_trace_as_reference: bool,
}

impl IsogenyGraphNodeFrobeniusData {
    /// Returns the graph node identifier.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the Frobenius trace package of this node representative.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns whether this node has the same order as the chosen reference node.
    pub fn same_curve_order_as_reference(&self) -> bool {
        self.same_curve_order_as_reference
    }

    /// Returns whether this node has the same Frobenius trace as the chosen reference node.
    pub fn same_trace_as_reference(&self) -> bool {
        self.same_trace_as_reference
    }

    /// Returns whether both comparisons against the reference node hold.
    ///
    /// Complexity: `Θ(1)`.
    pub fn holds(&self) -> bool {
        self.same_curve_order_as_reference && self.same_trace_as_reference
    }
}

/// Frobenius invariance report across the node representatives of an isogeny graph.
///
/// The current graph layer stores one representative curve at each node.
/// This report computes the Frobenius trace of each representative and compares
/// every node against one distinguished reference node.
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
    /// Returns the node chosen as the comparison reference.
    pub fn reference_node(&self) -> IsogenyGraphNodeId {
        self.reference_node
    }

    /// Returns the Frobenius trace package of the reference node.
    pub fn reference(&self) -> &FrobeniusTrace {
        &self.reference
    }

    /// Returns how many node representatives were checked.
    pub fn checked_nodes(&self) -> usize {
        self.checked_nodes
    }

    /// Returns how many stored graph edges accompany those node representatives.
    pub fn checked_edges(&self) -> usize {
        self.checked_edges
    }

    /// Returns the per-node Frobenius comparison data.
    pub fn node_reports(&self) -> &[IsogenyGraphNodeFrobeniusData] {
        &self.node_reports
    }

    /// Returns whether every checked node has the same order as the reference node.
    pub fn all_same_curve_order(&self) -> bool {
        self.all_same_curve_order
    }

    /// Returns whether every checked node has the same Frobenius trace as the reference node.
    pub fn all_same_trace(&self) -> bool {
        self.all_same_trace
    }

    /// Returns whether all checked node representatives agree with the reference.
    ///
    /// Complexity: `Θ(1)`.
    pub fn holds(&self) -> bool {
        self.all_same_curve_order && self.all_same_trace
    }
}

/// Exhaustively verifies the Frobenius relation for an explicit isogeny.
///
/// If the domain and codomain are curves over the same finite field `F_q`,
/// then one expects `#E(F_q) = #E'(F_q)` and equivalently the same Frobenius
/// trace `t_E = t_{E'}`.
///
/// The current implementation computes those traces by exhaustive point counts
/// on the two curves and compares the resulting invariants.
///
/// Complexity:
/// If `q = #F_q`, then this performs two exhaustive rational-point counts,
/// one on the domain and one on the codomain.
///
/// On the current `EnumerableCurveModel` path, one point count enumerates all
/// `x ∈ F_q` and tries to lift each `x` to zero, one, or two points, so the
/// mathematical cost is `Θ(q · C_lift) = `Θ(q²)`.
pub fn verify_isogeny_frobenius_relation<I, Domain, Codomain>(
    isogeny: &I,
) -> Result<IsogenyFrobeniusRelation, CurveError>
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
    let domain = isogeny.domain().frobenius_trace()?;
    let codomain = isogeny.codomain().frobenius_trace()?;

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
    Ok(IsogenyFrobeniusRelation {
        degree: isogeny.degree(),
        domain,
        codomain,
        same_curve_order,
        same_trace,
    })
}

/// Exhaustively verifies Frobenius invariance across one stored isogeny graph,
/// as expected for a graph built from curves that are mutually isogenous over
/// the same finite field.
///
/// Complexity:
/// If the graph has `V` nodes and the common base field has size `q`, this
/// performs `Θ(V)` exhaustive rational-point counts, one for each representative.
///
/// On the current `EnumerableCurveModel` path, one point count enumerates all
/// `x ∈ F_q` and tries to lift each `x` to zero, one, or two points, so the
/// mathematical cost is `Θ(Vq · C_lift) = `Θ(Vq²)`.
pub fn verify_isogeny_graph_frobenius_relation<C>(
    graph: &IsogenyGraph<C>,
) -> Result<IsogenyGraphFrobeniusReport, IsogenyGraphError>
where
    C: FrobeniusTraceCurveModel + GraphCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + core::hash::Hash + PartialEq,
    C::IsomorphismWitness: Clone + core::fmt::Debug,
{
    let reference_node = graph
        .nodes()
        .first()
        .ok_or(IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(0)))?;
    let reference = reference_node.representative().frobenius_trace()?;

    let mut node_reports = Vec::with_capacity(graph.node_count());
    let mut all_same_curve_order = true;
    let mut all_same_trace = true;

    for node in graph.nodes() {
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

        node_reports.push(IsogenyGraphNodeFrobeniusData {
            node_id: node.id(),
            frobenius_trace,
            same_curve_order_as_reference,
            same_trace_as_reference,
        });
    }

    Ok(IsogenyGraphFrobeniusReport {
        reference_node: reference_node.id(),
        reference,
        checked_nodes: graph.node_count(),
        checked_edges: graph.edge_count(),
        node_reports,
        all_same_curve_order,
        all_same_trace,
    })
}
