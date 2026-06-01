use core::marker::PhantomData;

use crate::elliptic_curves::{CurveError, CurveIsomorphism, CurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

use super::{Isogeny, IsogenyError};

type CurveBase<C> = <C as CurveModel>::BaseField;
type CurveElem<C> = <C as CurveModel>::Elem;
type CurvePoint<C> = <C as CurveModel>::Point;

pub trait CompositionBridge<Middle: CurveModel> {
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError>;

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError>;
}

impl<Middle: CurveModel + PartialEq> CompositionBridge<Middle> for () {
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError> {
        if source != target {
            return Err(IsogenyError::CompositionDomainCodomainMismatch);
        }

        Ok(())
    }

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError> {
        Ok(point)
    }
}

impl<Middle, Bridge> CompositionBridge<Middle> for Bridge
where
    Middle: CurveModel + Clone + PartialEq,
    Bridge: CurveIsomorphism<Domain = Middle, Codomain = Middle>,
{
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError> {
        if self.domain() != source || self.codomain() != target {
            return Err(IsogenyError::CompositionDomainCodomainMismatch);
        }

        Ok(())
    }

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError> {
        self.evaluate(&point)
            .map_err(|_| IsogenyError::ImagePointNotOnCodomain)
    }
}

/// Formal composition `second ∘ first` of two explicit isogenies.
///
/// This educational scaffold models the usual right-to-left composition:
///
/// - `first : E -> E'`
/// - `second : E' -> E''`
/// - `second ∘ first : E -> E''`
///
/// The field names intentionally use `first` and `second` instead of
/// `left` and `right`, since composition is read right-to-left and those names
/// are often easier to follow in educational code and prose.
///
/// At the current milestone this type validates that the middle curves agree
/// exactly at construction time and already supports composed point
/// evaluation. The remaining interface is left as explicit `todo!()`
/// placeholders for the next implementation step.
pub struct ComposedIsogeny<I, J, Domain, Middle, Codomain, Bridge = ()>
where
    Domain: FiniteGroupCurveModel,
    CurveBase<Domain>:
        EnumerableFiniteField<Elem = CurveElem<Domain>> + SqrtField<Elem = CurveElem<Domain>>,
    CurvePoint<Domain>: Clone + PartialEq,
    Middle: CurveModel,
    Codomain: CurveModel,
    CurvePoint<Codomain>: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    first: I,
    bridge: Bridge,
    second: J,
    kernel_points: Vec<Domain::Point>,
    marker: PhantomData<(Domain, Middle, Codomain)>,
}

impl<I, J, Domain, Middle, Codomain> ComposedIsogeny<I, J, Domain, Middle, Codomain, ()>
where
    Domain: FiniteGroupCurveModel,
    CurveBase<Domain>:
        EnumerableFiniteField<Elem = CurveElem<Domain>> + SqrtField<Elem = CurveElem<Domain>>,
    CurvePoint<Domain>: Clone + PartialEq,
    Middle: CurveModel + PartialEq,
    Codomain: CurveModel,
    CurvePoint<Codomain>: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
{
    /// Returns the first map in the composition.
    pub fn first(&self) -> &I {
        &self.first
    }

    /// Returns the second map in the composition.
    pub fn second(&self) -> &J {
        &self.second
    }

    /// Builds the formal composition `second ∘ first` under strict middle-curve equality.
    ///
    /// Since the public [`Isogeny`] trait exposes explicit domain and codomain
    /// curve values, this constructor first checks that the codomain carried by
    /// `first` is exactly the same curve value as the domain carried by
    /// `second`. If the middle curves differ, composition is rejected with
    /// [`IsogenyError::CompositionDomainCodomainMismatch`].
    ///
    /// This constructor is intentionally strict: it does not search for, infer,
    /// or apply any bridging isomorphism between the middle curves.
    pub fn new_strict(first: I, second: J) -> Result<Self, IsogenyError> {
        ().validate_bridge(first.codomain(), second.domain())?;
        let bridge = ();
        let kernel_points = compute_composed_kernel_points(&first, &bridge, &second)?;

        Ok(Self {
            first,
            bridge,
            second,
            kernel_points,
            marker: PhantomData,
        })
    }
}

impl<I, J, Domain, Middle, Codomain, Bridge> ComposedIsogeny<I, J, Domain, Middle, Codomain, Bridge>
where
    Domain: FiniteGroupCurveModel,
    CurveBase<Domain>:
        EnumerableFiniteField<Elem = CurveElem<Domain>> + SqrtField<Elem = CurveElem<Domain>>,
    CurvePoint<Domain>: Clone + PartialEq,
    Middle: CurveModel + Clone + PartialEq,
    Codomain: CurveModel,
    CurvePoint<Codomain>: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
    Bridge: CurveIsomorphism<Domain = Middle, Codomain = Middle>,
{
    /// Builds the bridged composition `second ∘ bridge ∘ first`.
    ///
    /// This constructor is for the common educational situation where
    /// `codomain(first)` and `domain(second)` are not literally the same curve
    /// equation, but are linked by an explicit bridge isomorphism
    ///
    /// `bridge : codomain(first) -> domain(second)`.
    ///
    /// The resulting composition therefore represents
    ///
    /// `second ∘ bridge ∘ first`.
    pub fn new_up_to_isomorphism(
        first: I,
        bridge: Bridge,
        second: J,
    ) -> Result<Self, IsogenyError> {
        bridge.validate_bridge(first.codomain(), second.domain())?;
        let kernel_points = compute_composed_kernel_points(&first, &bridge, &second)?;

        Ok(Self {
            first,
            bridge,
            second,
            kernel_points,
            marker: PhantomData,
        })
    }
}

impl<I, J, Domain, Middle, Codomain, Bridge> Isogeny<Domain, Codomain>
    for ComposedIsogeny<I, J, Domain, Middle, Codomain, Bridge>
where
    Domain: FiniteGroupCurveModel,
    CurveBase<Domain>:
        EnumerableFiniteField<Elem = CurveElem<Domain>> + SqrtField<Elem = CurveElem<Domain>>,
    CurvePoint<Domain>: Clone + PartialEq,
    Middle: CurveModel + Clone + PartialEq,
    Codomain: CurveModel,
    CurvePoint<Codomain>: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    /// Returns the domain of the composed map, namely the domain of `first`.
    fn domain(&self) -> &Domain {
        self.first.domain()
    }

    /// Returns the codomain of the composed map, namely the codomain of `second`.
    fn codomain(&self) -> &Codomain {
        self.second.codomain()
    }

    /// Returns the degree of the formal composition.
    ///
    /// `deg(second ∘ bridge ∘ first) = deg(first) * deg(second)`.
    ///
    /// Since [`ComposedIsogeny`] stores the two component maps explicitly, this
    /// method reads the two isogeny degrees directly and returns the plain
    /// product. The bridge is an isomorphism, so it contributes degree `1`.
    fn degree(&self) -> usize {
        self.first.degree() * self.second.degree()
    }

    /// Evaluates the composition at a point of the first domain.
    ///
    /// This follows the usual right-to-left rule
    ///
    /// `second ∘ bridge ∘ first (P) = second(bridge(first(P)))`.
    ///
    /// The current educational implementation also performs a few explicit
    /// confidence checks around that flow:
    ///
    /// - it rejects inputs that do not lie on the domain of `first`
    /// - it checks that the intermediate image returned by `first` really lies
    ///   on the codomain of `first`
    /// - it transports that image across the stored bridge, which is either
    ///   the strict identity bridge or an explicit curve isomorphism
    /// - it checks that the bridged image lies on the domain of `second`
    /// - it checks that the final image returned by `second` lies on the
    ///   declared final codomain
    fn evaluate(&self, point: &Domain::Point) -> Result<Codomain::Point, IsogenyError> {
        evaluate_composed_point(&self.first, &self.bridge, &self.second, point)
    }

    fn kernel_points(&self) -> &[Domain::Point] {
        &self.kernel_points
    }
}

fn compute_composed_kernel_points<Domain, Middle, Codomain, First, Second, Bridge>(
    first: &First,
    bridge: &Bridge,
    second: &Second,
) -> Result<Vec<Domain::Point>, IsogenyError>
where
    Domain: FiniteGroupCurveModel,
    CurveBase<Domain>:
        EnumerableFiniteField<Elem = CurveElem<Domain>> + SqrtField<Elem = CurveElem<Domain>>,
    CurvePoint<Domain>: Clone + PartialEq,
    Middle: CurveModel,
    Codomain: CurveModel,
    CurvePoint<Codomain>: PartialEq,
    First: Isogeny<Domain, Middle>,
    Second: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    let codomain_identity = second.codomain().identity();

    let kernel_points = first
        .domain()
        .points()
        .into_iter()
        .map(|point| {
            let image = evaluate_composed_point(first, bridge, second, &point)?;
            Ok((point, image == codomain_identity))
        })
        .collect::<Result<Vec<_>, IsogenyError>>()?
        .into_iter()
        .filter_map(|(point, maps_to_identity)| maps_to_identity.then_some(point))
        .collect::<Vec<_>>();

    Ok(kernel_points)
}

fn evaluate_composed_point<Domain, Middle, Codomain, First, Second, Bridge>(
    first: &First,
    bridge: &Bridge,
    second: &Second,
    point: &Domain::Point,
) -> Result<Codomain::Point, IsogenyError>
where
    Domain: CurveModel,
    Middle: CurveModel,
    Codomain: CurveModel,
    First: Isogeny<Domain, Middle>,
    Second: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    if !first.domain().contains(point) {
        return Err(CurveError::PointNotOnCurve.into());
    }

    let middle = first.evaluate(point)?;
    if !first.codomain().contains(&middle) {
        return Err(IsogenyError::ImagePointNotOnCodomain);
    }

    let bridged = bridge.transport(middle)?;
    if !second.domain().contains(&bridged) {
        return Err(IsogenyError::ImagePointNotOnCodomain);
    }

    let image = second.evaluate(&bridged)?;
    if !second.codomain().contains(&image) {
        return Err(IsogenyError::ImagePointNotOnCodomain);
    }

    Ok(image)
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, CurveError, CurveIsomorphism, CurveModel,
        EnumerableCurveModel, FiniteGroupCurveModel, ShortWeierstrassCurve,
        ShortWeierstrassIsomorphism,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::{
        ComposedIsogeny, Isogeny, IsogenyError, ScalarMultiplicationIsogeny, VeluIsogeny,
        VerifiableIsogeny, maps_equal_exhaustively,
    };

    type F41 = Fp<41>;
    type Curve = ShortWeierstrassCurve<F41>;
    type Velu = VeluIsogeny<Curve>;

    fn curve_a() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn first_generator(curve: &Curve) -> AffinePoint<F41> {
        curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("known two-torsion point should lie on the curve")
    }

    fn small_nontrivial_generator(curve: &Curve) -> AffinePoint<F41> {
        curve
            .point_orders()
            .into_iter()
            .filter(|(point, order)| *order > 1 && !curve.is_identity(point))
            .min_by_key(|(_, order)| *order)
            .map(|(point, _)| point)
            .expect("small sample curve should have a non-trivial point")
    }

    fn first_isogeny() -> Velu {
        let domain = curve_a();
        VeluIsogeny::from_generator(domain.clone(), first_generator(&domain))
            .expect("first sample Vélu isogeny should build")
    }

    fn second_isogeny(domain: &Curve) -> Velu {
        VeluIsogeny::from_generator(domain.clone(), small_nontrivial_generator(domain))
            .expect("second sample Vélu isogeny should build")
    }

    fn bridged_second_isogeny(first: &Velu) -> (ShortWeierstrassIsomorphism<F41>, Velu) {
        let bridge = ShortWeierstrassIsomorphism::new(first.codomain().clone(), F41::from_i64(3))
            .expect("sample bridge isomorphism should build");
        let generator = small_nontrivial_generator(first.codomain());
        let transported_generator = bridge
            .evaluate(&generator)
            .expect("bridge should transport the sample generator");
        let second = VeluIsogeny::from_generator(bridge.codomain().clone(), transported_generator)
            .expect("bridged second sample Vélu isogeny should build");

        (bridge, second)
    }

    fn first_non_kernel_point(isogeny: &Velu) -> AffinePoint<F41> {
        isogeny
            .domain()
            .points()
            .into_iter()
            .find(|point| !isogeny.kernel_points().contains(point))
            .expect("sample Vélu isogeny should have at least one point outside its kernel")
    }

    struct BrokenMiddleImageIsogeny {
        inner: Velu,
        broken_point: AffinePoint<F41>,
    }

    impl Isogeny<Curve, Curve> for BrokenMiddleImageIsogeny {
        fn domain(&self) -> &Curve {
            self.inner.domain()
        }

        fn codomain(&self) -> &Curve {
            self.inner.codomain()
        }

        fn degree(&self) -> usize {
            self.inner.degree()
        }

        fn evaluate(
            &self,
            point: &<Curve as CurveModel>::Point,
        ) -> Result<<Curve as CurveModel>::Point, IsogenyError> {
            if *point == self.broken_point {
                return Ok(AffinePoint::new(F41::from_i64(2), F41::from_i64(2)));
            }

            self.inner.evaluate(point)
        }

        fn kernel_points(&self) -> &[<Curve as CurveModel>::Point] {
            self.inner.kernel_points()
        }
    }

    #[test]
    fn new_strict_accepts_exact_matching_middle_curves() {
        let first = first_isogeny();
        let second = second_isogeny(first.codomain());

        let composed =
            ComposedIsogeny::new_strict(first, second).expect("composition should build");

        assert_eq!(composed.first().codomain(), composed.second().domain());
    }

    #[test]
    fn new_strict_rejects_mismatched_middle_curves() {
        let first = first_isogeny();
        let unrelated_second = second_isogeny(&curve_a());

        assert!(matches!(
            ComposedIsogeny::new_strict(first, unrelated_second),
            Err(IsogenyError::CompositionDomainCodomainMismatch)
        ));
    }

    #[test]
    fn degree_multiplies_the_component_degrees() {
        let first = first_isogeny();
        let second = second_isogeny(first.codomain());
        let expected_degree = first.degree() * second.degree();

        let composed =
            ComposedIsogeny::new_strict(first, second).expect("composition should build");

        assert_eq!(composed.degree(), expected_degree);
    }

    #[test]
    fn evaluate_applies_first_then_second_on_the_shared_middle_curve() {
        let first = first_isogeny();
        let second = second_isogeny(first.codomain());
        let point = first_non_kernel_point(&first);
        let expected = second
            .evaluate(
                &first
                    .evaluate(&point)
                    .expect("first Vélu isogeny should evaluate on sample point"),
            )
            .expect("second Vélu isogeny should evaluate on the transported point");

        let composed =
            ComposedIsogeny::new_strict(first, second).expect("composition should build");

        assert_eq!(
            composed
                .evaluate(&point)
                .expect("composition should evaluate"),
            expected
        );
    }

    #[test]
    fn evaluate_rejects_points_outside_the_first_domain_before_delegating() {
        let first = first_isogeny();
        let second = second_isogeny(first.codomain());
        let invalid = AffinePoint::new(F41::from_i64(2), F41::from_i64(2));

        let composed =
            ComposedIsogeny::new_strict(first, second).expect("composition should build");

        assert_eq!(
            composed.evaluate(&invalid),
            Err(IsogenyError::Curve(CurveError::PointNotOnCurve))
        );
    }

    #[test]
    fn constructor_rejects_intermediate_images_outside_the_shared_middle_curve() {
        let inner_first = first_isogeny();
        let point = first_non_kernel_point(&inner_first);
        let broken_first = BrokenMiddleImageIsogeny {
            inner: inner_first,
            broken_point: point.clone(),
        };
        let second = second_isogeny(broken_first.codomain());

        assert!(matches!(
            ComposedIsogeny::new_strict(broken_first, second),
            Err(IsogenyError::ImagePointNotOnCodomain)
        ));
    }

    #[test]
    fn new_up_to_isomorphism_accepts_a_valid_bridge_and_evaluates_as_psi_alpha_phi() {
        let first = first_isogeny();
        let point = first_non_kernel_point(&first);
        let middle_image = first
            .evaluate(&point)
            .expect("first Vélu isogeny should evaluate on sample point");
        let (bridge, second) = bridged_second_isogeny(&first);
        let expected = second
            .evaluate(
                &bridge
                    .evaluate(&middle_image)
                    .expect("bridge should transport the middle image"),
            )
            .expect("second Vélu isogeny should evaluate on the bridged point");

        let composed = ComposedIsogeny::new_up_to_isomorphism(first, bridge, second)
            .expect("bridged composition should build");

        assert_eq!(
            composed
                .evaluate(&point)
                .expect("bridged composition should evaluate"),
            expected
        );
    }

    #[test]
    fn new_up_to_isomorphism_rejects_a_bridge_with_the_wrong_domain() {
        let first = first_isogeny();
        let wrong_bridge = ShortWeierstrassIsomorphism::new(curve_a(), F41::from_i64(3))
            .expect("sample wrong bridge should still be a valid isomorphism");
        let second = VeluIsogeny::from_generator(
            wrong_bridge.codomain().clone(),
            small_nontrivial_generator(wrong_bridge.codomain()),
        )
        .expect("second Vélu isogeny on the wrong bridge codomain should build");

        assert!(matches!(
            ComposedIsogeny::new_up_to_isomorphism(first, wrong_bridge, second),
            Err(IsogenyError::CompositionDomainCodomainMismatch)
        ));
    }

    #[test]
    fn strict_composition_kernel_points_match_the_full_identity_fiber() {
        let first = first_isogeny();
        let second = second_isogeny(first.codomain());
        let composed =
            ComposedIsogeny::new_strict(first, second).expect("composition should build");

        assert_eq!(composed.verify_kernel_exactness(), Ok(()));
    }

    #[test]
    fn composition_with_identity_like_scalar_one_is_neutral() {
        let first = first_isogeny();
        let identity_like = ScalarMultiplicationIsogeny::new(first.codomain().clone(), 1)
            .expect("multiplication by one should build");
        let composed = ComposedIsogeny::new_strict(first.clone(), identity_like)
            .expect("composition with scalar-one should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve, Curve>(&composed, &first),
            Ok(true)
        );
    }
}
