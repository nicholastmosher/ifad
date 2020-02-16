use crate::{Aspect, AnnotationStatus, AnnotationRecord, GeneRecord};
use crate::index::AnnoIndex;

#[derive(Debug, Hash, Eq, PartialEq)]
pub struct Annotation {
    pub db: String,
    pub database_id: String,
    pub db_object_symbol: String,
    pub invert: bool,
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

impl Annotation {
    pub fn from_record(record: AnnotationRecord, experimental_evidence: &[&str]) -> Self {
        let mut gene_names = vec![record.unique_gene_name];
        gene_names.extend(record.alternative_gene_name.split('|').map(ToOwned::to_owned));

        let annotation_status =
            if &record.evidence_code == "ND" {
                AnnotationStatus::Unknown
            } else if experimental_evidence.contains(&&*record.evidence_code) {
                AnnotationStatus::KnownExperimental
            } else {
                AnnotationStatus::KnownOther
            };

        Annotation {
            db: record.db,
            database_id: record.database_id,
            db_object_symbol: record.db_object_symbol,
            invert: record.invert.to_lowercase() == "not",
            go_term: record.go_term,
            reference: record.reference,
            evidence_code: record.evidence_code,
            additional_evidence: record.additional_evidence,
            aspect: record.aspect,
            annotation_status,
            gene_names,
            gene_product_type: record.gene_product_type,
            taxon: record.taxon,
            date: record.date,
            assigned_by: record.assigned_by,
            annotation_extension: record.annotation_extension,
            gene_product_form_id: record.gene_product_form_id,
        }
    }

    pub fn gene_in(&self, index: &AnnoIndex) -> Option<&str> {
        self.gene_names.iter()
            .find(|name| index.contains_key(*name))
            .map(String::as_str)
    }
}

#[derive(Debug, Hash, Eq, PartialEq)]
pub struct Gene {
    pub gene_id: String,
    pub gene_product_type: String,
}

impl From<GeneRecord> for Gene {
    fn from(record: GeneRecord) -> Self {
        Gene {
            gene_id: record.gene_id,
            gene_product_type: record.gene_product_type,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convert_annotation() {
        let record = AnnotationRecord {
            db: "TAIR".to_string(),
            database_id: "locus:2031476".to_string(),
            db_object_symbol: "ENO1".to_string(),
            invert: "".to_string(),
            go_term: "GO:0000015".to_string(),
            reference: "TAIR:AnalysisReference:501756966".to_string(),
            evidence_code: "IEA".to_string(),
            additional_evidence: "InterPro:IPR000941".to_string(),
            aspect: Aspect::CellularComponent,
            unique_gene_name: "AT1G74030".to_string(),
            alternative_gene_name: "AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10".to_string(),
            gene_product_type: "protein".to_string(),
            taxon: "taxon:3702".to_string(),
            date: "20190907".to_string(),
            assigned_by: "InterPro".to_string(),
            annotation_extension: "".to_string(),
            gene_product_form_id: "TAIR:locus:2031476".to_string(),
        };
        let annotation = Annotation::from_record(record, &["IEA"]);
        let expected_annotation = Annotation {
            db: "TAIR".to_string(),
            database_id: "locus:2031476".to_string(),
            db_object_symbol: "ENO1".to_string(),
            invert: false,
            go_term: "GO:0000015".to_string(),
            reference: "TAIR:AnalysisReference:501756966".to_string(),
            evidence_code: "IEA".to_string(),
            additional_evidence: "InterPro:IPR000941".to_string(),
            aspect: Aspect::CellularComponent,
            annotation_status: AnnotationStatus::KnownExperimental,
            gene_names: vec![
                "AT1G74030".to_string(),
                "AT1G74030".to_string(),
                "ENO1".to_string(),
                "enolase 1".to_string(),
                "F2P9.10".to_string(),
                "F2P9_10".to_string(),
            ],
            gene_product_type: "protein".to_string(),
            taxon: "taxon:3702".to_string(),
            date: "20190907".to_string(),
            assigned_by: "InterPro".to_string(),
            annotation_extension: "".to_string(),
            gene_product_form_id: "TAIR:locus:2031476".to_string(),
        };
        assert_eq!(annotation, expected_annotation);
    }
}
