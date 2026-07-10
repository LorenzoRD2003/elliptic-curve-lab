/// Purely graph-theoretic reason why a deterministic crater walk stopped.
///
/// These variants describe only the stored crater evidence and the deterministic
/// edge-ordering rule. They do not certify an arithmetic orientation, and they
/// should not be compared with the order of an ideal class without additional
/// orientation data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CraterWalkTermination {
    /// The path returned to the requested starting node.
    ClosedCycle,
    /// The requested starting node was not among the certified crater nodes.
    StartOutsideCrater,
    /// The current crater node had no certified outgoing internal crater edge.
    NoCertifiedOutgoingEdge,
    /// The walk encountered a previously visited node different from `start`.
    RepeatedNonStartNode,
}
