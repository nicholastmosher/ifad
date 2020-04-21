use crate::{Aspect, AnnotationStatus, AnnotationRecord, GeneRecord};
use crate::index::{GeneKey, AnnoIndex};

#[derive(Debug, Hash, Eq, PartialEq, Clone)]
pub struct Annotation {
    pub record: AnnotationRecord,
    pub invert: bool,
    pub aspect: Aspect,
    pub annotation_status: AnnotationStatus,
}

impl Annotation {
    pub fn from_record(record: AnnotationRecord, experimental_evidence: &[&str]) -> Annotation {
        let annotation_status =
            if &record.evidence_code == "ND" {
                AnnotationStatus::Unknown
            } else if experimental_evidence.contains(&&*record.evidence_code) {
                AnnotationStatus::KnownExperimental
            } else {
                AnnotationStatus::KnownOther
            };

        Annotation {
            invert: record.invert.eq_ignore_ascii_case("not"),
            aspect: record.aspect,
            annotation_status,
            record
        }
    }

    pub fn gene_names(&self) -> impl Iterator<Item=&str> {
        std::iter::once(&*self.record.unique_gene_name)
            .chain(self.record.alternative_gene_name.split('|'))
    }

    pub fn gene_in(&self, index: &AnnoIndex) -> Option<GeneKey> {
        self.gene_names()
            .find(|name| index.contains_key(&((**name).to_string())))
            .and_then(|name| index.get(name).map(|(gene, _)| *gene))
    }
}

#[derive(Debug, Hash, Eq, PartialEq, Clone)]
pub struct Gene {
    pub record: GeneRecord,
}

impl Gene {
    pub fn from_record(record: GeneRecord) -> Gene {
        Gene { record }
    }

    #[inline(always)]
    pub fn gene_id(&self) -> &str {
        &self.record.gene_id
    }

    #[inline(always)]
    pub fn gene_product_type(&self) -> &str {
        &self.record.gene_product_type
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
        let annotation = Annotation::from_record(record.clone(), &["IEA"]);
        let expected_annotation = Annotation {
            invert: false,
            aspect: Aspect::CellularComponent,
            annotation_status: AnnotationStatus::KnownExperimental,
            record,
        };
        assert_eq!(annotation, expected_annotation);
    }
}
