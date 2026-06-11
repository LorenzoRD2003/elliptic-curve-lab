mod absolute;
mod relative;
mod shared;
mod r#trait;
mod verschiebung;

pub use absolute::AbsoluteFrobeniusIsogeny;
pub use relative::RelativeFrobeniusIsogeny;
pub use r#trait::FrobeniusLikeIsogeny;
pub use verschiebung::{VerschiebungCertificate, VerschiebungIsogeny};

#[cfg(test)]
mod tests;
