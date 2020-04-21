use futures::FutureExt;
use actix_web::{HttpServer, App, web, Responder};
use ifad::{MetadataReader, Gene, Annotation, Index, Query, StreamingGafExporter};
use std::sync::Arc;
use arc_swap::ArcSwap;
use actix_web::dev::HttpResponseBuilder;
use actix_web::http::StatusCode;
use actix_web::error::PayloadError;
use std::io::BufReader;

struct Config {
    genes_file: String,
    annotations_file: String,
}

impl Config {
    fn from_env() -> Option<Config> {
        let genes_file = std::env::var("GENES_FILE").ok()?;
        let annotations_file = std::env::var("ANNOTATIONS_FILE").ok()?;
        Some(Config { genes_file, annotations_file })
    }
}

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
    // let gene_metadata = gene_reader.metadata().expect("should capture gene metadata").to_string();
    // let gene_headers = gene_reader.header().expect("should get gene headers").to_string();
    let genes: Vec<Gene> = gene_records.into_iter()
        .map(|record| Gene::from_record(record))
        .collect();

    let mut annos_file = std::fs::File::open(config.annotations_file)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    // let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata").to_string();
    // let anno_headers = anno_reader.header().expect("should capture annotation header").to_string();
    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_records.into_iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();

    let index = Index::new(genes, annotations);
    let swap = ArcSwap::new(Arc::new(index));
    actix::System::new("ifad").block_on(server(swap).map(|_| ()));
    Ok(())
}

async fn index(data: web::Data<ArcSwap<Index>>) -> impl Responder {
    use futures::StreamExt;

    let index = data.get_ref().load_full();
    let query_result = Query::All.execute(index);
    let stream = StreamingGafExporter::new(
        "anno metadata".to_string(),
        "anno_header".to_string(),
        query_result.iter_annotations().map(|anno| anno.record),
    ).map(|result| result.map_err(|_| PayloadError::EncodingCorrupted));

    HttpResponseBuilder::new(StatusCode::OK)
        .streaming(stream)
}

async fn server(data: ArcSwap<Index>) -> std::io::Result<()> {
    HttpServer::new(move || {
        App::new()
            .data(data.clone())
            .route("/", web::get().to(index))
    })
    .bind("127.0.0.1:8000")?
    .run()
    .await
}
