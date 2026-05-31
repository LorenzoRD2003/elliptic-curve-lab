/// Common interface for curve models.
pub trait CurveEquation<Point> {
    /// Returns whether the given point satisfies the curve equation.
    fn is_on_curve(&self, point: &Point) -> bool;
}
