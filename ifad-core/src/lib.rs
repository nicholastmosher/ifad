#![deny(warnings)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

mod ingest;
mod models;
mod index;
mod queries;

pub use ingest::{AnnotationRecord, GeneRecord, MetadataReader};
pub use models::{Annotation, Gene};
pub use index::Index;

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum Aspect {
    #[serde(rename = "F")]
    MolecularFunction,
    #[serde(rename = "P")]
    BiologicalProcess,
    #[serde(rename = "C")]
    CellularComponent,
}

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum AnnotationStatus {
    KnownExperimental,
    KnownOther,
    Unknown,
    Unannotated,
}
