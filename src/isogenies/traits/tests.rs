use crate::elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveError, CurveModel, EnumerableCurveModel,
    ShortWeierstrassCurve,
};
use crate::fields::{Field, Fp};
use crate::isogenies::{Isogeny, IsogenyError, VeluIsogeny, VerifiableIsogeny};

type F41 = Fp<41>;
type Curve = ShortWeierstrassCurve<F41>;
type Point = AffinePoint<F41>;

fn f41_curve() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn sample_point(curve: &Curve) -> Point {
    curve
        .point(F41::from_i64(3), F41::from_i64(6))
        .expect("point should lie on the curve")
}

fn explicit_identity_kernel(curve: &Curve) -> Vec<Point> {
    vec![curve.identity()]
}

#[derive(Clone)]
struct TableIsogeny {
    domain: Curve,
    codomain: Curve,
    kernel_points: Vec<Point>,
    images: Vec<(Point, Point)>,
}

impl TableIsogeny {
    fn new(
        domain: Curve,
        codomain: Curve,
        kernel_points: Vec<Point>,
        image_fn: impl Fn(&Point) -> Point,
    ) -> Self {
        let images = domain
            .points()
            .into_iter()
            .map(|point| {
                let image = image_fn(&point);
                (point, image)
            })
            .collect();

        Self {
            domain,
            codomain,
            kernel_points,
            images,
        }
    }
}

impl Isogeny<Curve, Curve> for TableIsogeny {
    fn domain(&self) -> &Curve {
        &self.domain
    }

    fn codomain(&self) -> &Curve {
        &self.codomain
    }

    fn degree(&self) -> usize {
        self.kernel_points.len()
    }

    fn evaluate(&self, point: &Point) -> Result<Point, IsogenyError> {
        if !self.domain.contains(point) {
            return Err(CurveError::PointNotOnCurve.into());
        }

        self.images
            .iter()
            .find_map(|(source, image)| (source == point).then_some(image.clone()))
            .ok_or_else(|| CurveError::PointNotOnCurve.into())
    }

    fn kernel_points(&self) -> &[Point] {
        &self.kernel_points
    }
}

#[test]
fn velu_isogeny_passes_all_exhaustive_verifiers_on_the_small_f41_example() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny = VeluIsogeny::from_generator(domain, generator).expect("isogeny should build");

    assert_eq!(isogeny.verify_maps_domain_to_codomain(), Ok(()));
    assert_eq!(isogeny.verify_maps_kernel_to_identity(), Ok(()));
    assert_eq!(isogeny.verify_homomorphism(), Ok(()));
    assert_eq!(isogeny.verify_kernel_exactness(), Ok(()));
}

#[test]
fn verify_maps_domain_to_codomain_detects_an_image_outside_the_declared_codomain() {
    let curve = f41_curve();
    let bad_point = sample_point(&curve);
    let off_curve_image = AffinePoint::new(F41::from_i64(2), F41::from_i64(2));
    let isogeny = TableIsogeny::new(
        curve.clone(),
        curve.clone(),
        explicit_identity_kernel(&curve),
        |point| {
            if *point == bad_point {
                off_curve_image.clone()
            } else {
                point.clone()
            }
        },
    );

    assert_eq!(
        isogeny.verify_maps_domain_to_codomain(),
        Err(IsogenyError::ImagePointNotOnCodomain)
    );
}

#[test]
fn verify_maps_kernel_to_identity_detects_a_declared_kernel_point_with_nonzero_image() {
    let curve = f41_curve();
    let point = sample_point(&curve);
    let isogeny = TableIsogeny::new(
        curve.clone(),
        curve.clone(),
        vec![curve.identity(), point.clone()],
        |candidate| candidate.clone(),
    );

    assert_eq!(
        isogeny.verify_maps_kernel_to_identity(),
        Err(IsogenyError::KernelPointDoesNotMapToIdentity)
    );
}

#[test]
fn verify_homomorphism_detects_a_non_additive_mapping() {
    let curve = f41_curve();
    let bad_point = sample_point(&curve);
    let isogeny = TableIsogeny::new(
        curve.clone(),
        curve.clone(),
        vec![curve.identity(), bad_point.clone()],
        |point| {
            if *point == bad_point {
                curve.identity()
            } else {
                point.clone()
            }
        },
    );

    assert_eq!(
        isogeny.verify_homomorphism(),
        Err(IsogenyError::HomomorphismViolation)
    );
}

#[test]
fn verify_kernel_exactness_detects_when_the_declared_kernel_is_too_small() {
    let curve = f41_curve();
    let witness = sample_point(&curve);
    let isogeny = TableIsogeny::new(
        curve.clone(),
        curve.clone(),
        explicit_identity_kernel(&curve),
        |_| curve.identity(),
    );

    assert_eq!(
        isogeny
            .evaluate(&witness)
            .expect("constant map should evaluate"),
        curve.identity()
    );
    assert_eq!(
        isogeny.verify_kernel_exactness(),
        Err(IsogenyError::KernelMismatch)
    );
}
