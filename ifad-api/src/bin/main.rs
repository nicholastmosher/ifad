use std::sync::Arc;
use std::io::BufReader;
use futures::FutureExt;
use arc_swap::ArcSwap;
use ifad::{MetadataReader, Gene, Annotation, Index};

use ifad_api::app::{Config, AppData, server};

fn main() {
    match run() {
        Ok(()) => return,
        Err(e) => eprintln!("{}", e),
    }
}

fn run() -> Result<(), String> {
    dotenv::dotenv().map_err(|e| format!("failed to read .env: {:?}", e))?;
    let config = Config::from_env()
        .ok_or("failed to read Config from environment")?;

    let mut genes_file = std::fs::File::open(config.genes_file)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut gene_reader = MetadataReader::new(BufReader::new(&mut genes_file));
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}", e))?;
    let gene_metadata = gene_reader.metadata().expect("should capture gene metadata").to_string();
    let gene_headers = gene_reader.header().expect("should get gene headers").to_string();
    let genes: Vec<Gene> = gene_records.into_iter()
        .map(|record| Gene::from_record(record))
        .collect();

    let mut annos_file = std::fs::File::open(config.annotations_file)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata").to_string();
    let anno_headers = anno_reader.header().expect("should capture annotation header").to_string();
    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_records.into_iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();

    let index = Arc::new(Index::new(genes, annotations));
    let appdata = AppData {
        index,
        gene_metadata,
        gene_headers,
        anno_metadata,
        anno_headers,
    };
    let swap = ArcSwap::new(Arc::new(appdata));
    actix::System::new("ifad").block_on(server(swap).map(|_| ()));
    Ok(())
}
