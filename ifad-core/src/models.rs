use crate::{Aspect, AnnotationStatus};

#[derive(Debug, Hash, Eq, PartialEq)]
pub struct Annotation {
    pub db: String,
    pub database_id: String,
    pub db_object_symbol: String,
    pub invert: String,
    pub go_term: String,
    pub reference: String,
    pub evidence_code: String,
    pub additional_evidence: String,
    pub aspect: Aspect,
    pub annotation_status: AnnotationStatus,
    pub gene_names: Vec<String>,
    pub gene_product_type: String,
    pub taxon: String,
    pub date: String,
    pub assigned_by: String,
    pub annotation_extension: String,
    pub gene_product_form_id: String,
}

#[derive(Debug)]
pub struct Gene {
    pub gene_id: String,
    pub gene_product_type: String,
}
