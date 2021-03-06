use std::io::Write;
use serde::Serialize;

pub struct GafExporter<I: Iterator> {
    metadata: String,
    header: String,
    record_iter: I,
}

impl<T, I: Iterator<Item=T>> GafExporter<I>
    where T: Serialize
{
    pub fn new(
        metadata: String,
        header: String,
        record_iter: I,
    ) -> GafExporter<I> {
        GafExporter { metadata, header, record_iter }
    }

    pub fn write_all<W: Write>(&mut self, mut writer: W) -> std::io::Result<()> {
        write!(&mut writer, "{}", self.metadata)?;
        write!(&mut writer, "{}", self.header)?;

        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(writer);
        for record in &mut self.record_iter {
            csv_writer.serialize(record)?;
        }
        csv_writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::{AnnotationRecord, MetadataReader, GeneRecord};

    #[test]
    fn test_export_annotations() {
        let annotations_file = r"!gaf-version: 2.1
!
!Generated by GO Central
!
!Date Generated by GOC: 2019-10-07
DB	DB Object ID	DB Object Symbol	Qualifier	GO ID	DB:Reference (JDB:Reference)	Evidence Code	With (or) From	Aspect	DB Object Name	DB Object Type	Taxon	Date	Assigned By	Annotation Extension	Gene Product Form ID
TAIR	locus:2031476	ENO1		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT1G74030	AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10	protein	taxon:3702	20190907	InterPro		TAIR:locus:2031476
TAIR	locus:2043067	ENOC		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G29560	AT2G29560|ENOC|ENO3|cytosolic enolase|enolase 3|F16P2.6|F16P2_6	protein	taxon:3702	20190408	InterPro		TAIR:locus:2043067
TAIR	locus:2044851	LOS2		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G36530	AT2G36530|LOS2|ENO2|LOW EXPRESSION OF OSMOTICALLY RESPONSIVE GENES 2|enolase 2|F1O11.16|F1O11_16	protein	taxon:3702	20190408	InterPro		TAIR:locus:2044851
TAIR	locus:2032970	AT1G25260		GO:0000027	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR033867	P	AT1G25260	AT1G25260|F4F7.35|F4F7_35	protein	taxon:3702	20190404	InterPro		TAIR:locus:2032970
";
        let mut reader = MetadataReader::new(Cursor::new(&annotations_file));
        let annotations_records = AnnotationRecord::parse_from(&mut reader).expect("should parse annotations");
        let metadata = reader.metadata().expect("should get metadata");
        let header = reader.header().expect("should get header");

        // Test export
        let mut exporter = GafExporter::new(
            metadata.to_string(),
            header.to_string(),
            annotations_records.iter());

        let mut output = Vec::new();
        exporter.write_all(Cursor::new(&mut output)).unwrap();
        let output_string = String::from_utf8(output).unwrap();
        assert_eq!(&annotations_file, &output_string);
    }

    #[test]
    fn test_export_genes() {
        let genes_file = r"!Gene list based on the Araport11 genome release
!Release date of annotation: June 2016
!Annotated by: Araport at J. Craig Ventner Institute
name	gene_model_type
AT1G01010	protein_coding
AT1G01020	protein_coding
AT1G01030	protein_coding
AT1G01040	protein_coding
AT1G01046	miRNA_primary_transcript
";
        let mut reader = MetadataReader::new(Cursor::new(&genes_file));
        let genes_records = GeneRecord::parse_from(&mut reader).expect("should parse genes");
        let metadata = reader.metadata().expect("should get metadata");
        let header = reader.header().expect("should get header");

        // Test export
        let mut exporter = GafExporter::new(
            metadata.to_string(),
            header.to_string(),
            genes_records.iter());

        let mut output = Vec::new();
        exporter.write_all(Cursor::new(&mut output)).unwrap();
        let output_string = String::from_utf8(output).unwrap();
        assert_eq!(&genes_file, &output_string);
    }
}
