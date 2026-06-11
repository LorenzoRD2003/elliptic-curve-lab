mod absolute;
mod relative;
mod shared;
mod r#trait;

pub use absolute::AbsoluteFrobeniusIsogeny;
pub use relative::RelativeFrobeniusIsogeny;
pub use r#trait::FrobeniusLikeIsogeny;

#[cfg(test)]
mod tests;
