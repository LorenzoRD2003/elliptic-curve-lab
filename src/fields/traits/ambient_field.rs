/// Runtime-ambient field interface.
///
/// This trait is intentionally parallel to [`crate::fields::traits::Field`],
/// not a replacement for it.
///
/// The existing [`crate::fields::traits::Field`] trait is a static family
/// interface: its identities, arithmetic, characteristic, and small-integer
/// embedding are all determined by the implementing type alone.
///
/// Some mathematically natural fields in this repository depend on additional
/// runtime context. An example is the function field `F(E)` of one concrete
/// short-Weierstrass curve `E`, where arithmetic depends on the chosen cubic
/// relation `y^2 = x^3 + ax + b`.
///
/// For those ambient, runtime-dependent situations, this trait keeps the
/// algebraic operations on the family object itself:
///
/// - identities come from `&self`
/// - arithmetic can validate that two elements are in the same ambient family
/// - division can fail through the family's own error type
pub trait AmbientField {
    /// Stored element type of the ambient field.
    type Elem;

    /// Recoverable failure surface for ambient arithmetic.
    type Error;

    /// Returns the additive identity.
    fn zero(&self) -> Self::Elem;

    /// Returns the multiplicative identity.
    fn one(&self) -> Self::Elem;

    /// Returns whether two elements are equal in the ambient field.
    fn eq(&self, left: &Self::Elem, right: &Self::Elem) -> bool;

    /// Returns whether the given element is zero.
    fn is_zero(&self, value: &Self::Elem) -> bool {
        self.eq(value, &self.zero())
    }

    /// Returns whether the given element is one.
    fn is_one(&self, value: &Self::Elem) -> bool {
        self.eq(value, &self.one())
    }

    /// Adds two elements in the ambient field.
    fn add(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Returns the additive inverse of one element.
    fn neg(&self, value: &Self::Elem) -> Self::Elem;

    /// Subtracts two elements in the ambient field.
    fn sub(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        self.add(left, &self.neg(right))
    }

    /// Multiplies two elements in the ambient field.
    fn mul(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Returns the multiplicative inverse when it exists.
    fn inverse(&self, value: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Divides `left` by `right` when `right` is invertible.
    fn div(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        self.mul(left, &self.inverse(right)?)
    }
}
