#![deny(warnings)]
#![allow(dead_code)]

#[cfg(test)]
#[macro_use]
extern crate lazy_static;

use std::convert::TryFrom;
use serde::{Deserialize, Serialize};

mod ingest;
mod models;
mod index;
mod queries;
mod export;

pub use ingest::{AnnotationRecord, GeneRecord, MetadataReader};
pub use models::{Annotation, Gene};
pub use index::Index;
pub use queries::{Segment, Query, QueryResult};
pub use export::GafExporter;

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum Aspect {
    #[serde(rename = "F")]
    MolecularFunction,
    #[serde(rename = "P")]
    BiologicalProcess,
    #[serde(rename = "C")]
    CellularComponent,
}

impl TryFrom<&str> for Aspect {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let aspect = match value {
            "F" => Aspect::MolecularFunction,
            "P" => Aspect::MolecularFunction,
            "C" => Aspect::MolecularFunction,
            _ => return Err(()),
        };
        Ok(aspect)
    }
}

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum AnnotationStatus {
    KnownExperimental,
    KnownOther,
    Unknown,
    Unannotated,
}

impl TryFrom<&str> for AnnotationStatus {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let status = match value {
            "EXP" => AnnotationStatus::KnownExperimental,
            "OTHER" => AnnotationStatus::KnownOther,
            "UNKNOWN" => AnnotationStatus::Unknown,
            "UNANNOTATED" => AnnotationStatus::Unannotated,
            _ => return Err(()),
        };
        Ok(status)
    }
}
