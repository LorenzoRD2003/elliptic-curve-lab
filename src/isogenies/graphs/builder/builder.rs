use crate::isogenies::graphs::GraphCurveModel;

/// Small BFS-style builder for educational `ℓ`-isogeny graphs.
#[derive(Clone, Debug)]
pub struct IsogenyGraphBuilder<C: GraphCurveModel> {
    pub(crate) start_curve: C,
    pub(crate) ell: usize,
    pub(crate) max_depth: usize,
    pub(crate) deduplicate_by_base_field_isomorphism: bool,
}

impl<C: GraphCurveModel> IsogenyGraphBuilder<C> {
    /// Starts a graph build from one chosen representative and one prime degree `ℓ`.
    pub fn new(start_curve: C, degree: usize) -> Self {
        Self {
            start_curve,
            ell: degree,
            max_depth: 1,
            deduplicate_by_base_field_isomorphism: true,
        }
    }

    /// Sets the maximum BFS depth measured in edge traversals from the start node.
    pub fn max_depth(mut self, max_depth: usize) -> Self {
        self.max_depth = max_depth;
        self
    }

    /// Enables or disables base-field-isomorphism deduplication for newly
    /// discovered codomain curves.
    ///
    /// Exact representative equality is still deduplicated even when this flag
    /// is `false`.
    pub fn deduplicate_by_base_field_isomorphism(mut self, yes: bool) -> Self {
        self.deduplicate_by_base_field_isomorphism = yes;
        self
    }
}
