use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::{kernel::KernelDescription, traits::DegreeFactorizedIsogeny};
use num_bigint::BigUint;

/// High-level classification of the duality story represented by one report.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DualityKind {
    /// Classical duality for a separable small-field Vélu isogeny.
    SeparableClassical,
    /// Frobenius/Verschiebung duality in characteristic `p`.
    FrobeniusVerschiebung,
    /// Placeholder bucket for future partially modeled or mixed duality stories.
    MixedOrPartial,
}

/// Compact summary of separable/inseparable degree data.
///
/// When both entries are present, the intended interpretation is
///
/// `total degree = separable_degree * inseparable_degree`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DegreeFactorizationSummary {
    separable_degree: Option<BigUint>,
    inseparable_degree: Option<BigUint>,
}

/// Compact summary of the currently available kernel description.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KernelDescriptionSummary {
    pub(crate) total_degree: Option<usize>,
    pub(crate) reduced_degree: Option<usize>,
    pub(crate) infinitesimal_degree: Option<usize>,
    pub(crate) is_fully_reduced: bool,
    pub(crate) rational_point_count: Option<usize>,
    pub(crate) short_label: String,
}

impl DegreeFactorizationSummary {
    /// Builds a degree summary with no finer separable/inseparable data.
    pub fn unknown() -> Self {
        Self {
            separable_degree: None,
            inseparable_degree: None,
        }
    }

    /// Builds a degree summary from explicit optional factors.
    pub fn new(separable_degree: Option<BigUint>, inseparable_degree: Option<BigUint>) -> Self {
        Self {
            separable_degree,
            inseparable_degree,
        }
    }

    /// Builds a degree summary from an isogeny that already models the
    /// separable/inseparable factorization explicitly.
    pub fn from_degree_factorized_isogeny<I, D: CurveModel, C: CurveModel>(isogeny: &I) -> Self
    where
        I: DegreeFactorizedIsogeny<D, C>,
    {
        Self::new(
            Some(isogeny.separable_degree()),
            Some(isogeny.inseparable_degree()),
        )
    }

    /// Returns the separable-degree factor when known.
    pub fn separable_degree(&self) -> Option<&BigUint> {
        self.separable_degree.as_ref()
    }

    /// Returns the inseparable-degree factor when known.
    pub fn inseparable_degree(&self) -> Option<&BigUint> {
        self.inseparable_degree.as_ref()
    }
}

impl KernelDescriptionSummary {
    /// Builds a compact summary from a full kernel description.
    pub fn from_kernel_description<C: CurveModel>(description: &KernelDescription<C>) -> Self {
        let short_label = match description {
            KernelDescription::Reduced(_) => "reduced kernel".to_string(),
            KernelDescription::NonReduced(nonreduced) => {
                format!("nonreduced kernel ({})", nonreduced.label())
            }
            KernelDescription::Mixed(mixed) => match mixed.label() {
                Some(label) => format!("mixed kernel ({label})"),
                None => "mixed kernel".to_string(),
            },
            KernelDescription::Unknown => "unknown kernel description".to_string(),
        };

        Self {
            total_degree: description.degree(),
            reduced_degree: description.reduced_degree(),
            infinitesimal_degree: description.infinitesimal_degree(),
            is_fully_reduced: description.is_fully_reduced(),
            rational_point_count: description.rational_points().map(<[C::Point]>::len),
            short_label,
        }
    }

    /// Returns the currently known total kernel degree.
    pub fn total_degree(&self) -> Option<usize> {
        self.total_degree
    }

    /// Returns the currently known reduced kernel degree.
    pub fn reduced_degree(&self) -> Option<usize> {
        self.reduced_degree
    }

    /// Returns the currently known infinitesimal kernel degree.
    pub fn infinitesimal_degree(&self) -> Option<usize> {
        self.infinitesimal_degree
    }

    /// Returns whether the current kernel description is fully reduced.
    pub fn is_fully_reduced(&self) -> bool {
        self.is_fully_reduced
    }

    /// Returns how many explicit rational points are currently visible in the
    /// kernel description, when that notion makes sense for the current case.
    pub fn rational_point_count(&self) -> Option<usize> {
        self.rational_point_count
    }

    /// Returns the short human-readable kernel label used by reports and
    /// visualization helpers.
    pub fn short_label(&self) -> &str {
        &self.short_label
    }
}
