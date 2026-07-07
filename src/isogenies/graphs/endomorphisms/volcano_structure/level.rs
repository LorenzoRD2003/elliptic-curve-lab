use crate::isogenies::graphs::IsogenyGraphNodeId;

/// Role of a certified vertex or level in the observed `ℓ`-volcano structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VolcanoStructureRole {
    /// The single level of a depth-zero volcano, simultaneously `V_0` and `V_d`.
    SurfaceAndFloor,
    /// Level `0`, the surface `V_0`; this is also the crater in volcano language.
    Surface,
    /// An intermediate level strictly between the surface and the floor.
    Middle,
    /// The deepest certified level `V_d`.
    Floor,
}

impl VolcanoStructureRole {
    pub(crate) fn for_level(level: usize, certified_depth: usize) -> Self {
        if certified_depth == 0 {
            Self::SurfaceAndFloor
        } else if level == 0 {
            Self::Surface
        } else if level == certified_depth {
            Self::Floor
        } else {
            Self::Middle
        }
    }
}

/// One certified level in a [`super::VolcanoStructureReport`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VolcanoStructureLevelReport {
    level: usize,
    role: VolcanoStructureRole,
    nodes: Vec<IsogenyGraphNodeId>,
}

impl VolcanoStructureLevelReport {
    pub(crate) fn new(
        level: usize,
        role: VolcanoStructureRole,
        nodes: Vec<IsogenyGraphNodeId>,
    ) -> Self {
        Self { level, role, nodes }
    }

    /// Returns the level index, with `0` representing the surface `V₀`.
    pub fn level(&self) -> usize {
        self.level
    }

    /// Returns this level's structural role.
    pub fn role(&self) -> VolcanoStructureRole {
        self.role
    }

    /// Returns the stored nodes certified to lie on this level.
    pub fn nodes(&self) -> &[IsogenyGraphNodeId] {
        &self.nodes
    }
}
