use clap::{App, Arg, ArgMatches};
use std::io::BufReader;
use ifad::{MetadataReader, Annotation, Gene, Index, Segment, GafExporter, Query};
use std::convert::TryFrom;

fn app<'a, 'b>() -> clap::App<'a, 'b> {

    let segment_validator = |segment: String| {
        let splits: Vec<&str> = segment.split(",").collect();
        if splits.len() != 2 { return Err(format!("segments must be written as ASPECT,STATUS")); }
        if let &[aspect, status] = &splits[..] {
            if !&["F", "C", "P"].contains(&aspect) {
                return Err(format!("aspect must be one of F, C, or P"));
            }
            if !&["EXP", "OTHER", "UNKNOWN", "UNANNOTATED"].contains(&status) {
                return Err(format!("status must be one of EXP, OTHER, UNKNOWN, or UNANNOTATED"));
            }
            return Ok(());
        }
        unreachable!()
    };

    App::new("ifad")
        .arg(Arg::with_name("genes")
            .long("--genes")
            .required(true)
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations")
            .long("--annotations")
            .required(true)
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("genes_out")
            .long("--genes-out")
            .required(true)
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations_out")
            .long("--annotations-out")
            .required(true)
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("query")
            .long("--query")
            .possible_values(&["union", "intersection"])
            .default_value("union")
            .require_equals(true))
        .arg(Arg::with_name("segment")
            .multiple(true)
            .long("--segment")
            .required_unless("all")
            .require_equals(true)
            .takes_value(true)
            .validator(segment_validator))
}

fn main() {
    let matches = app().get_matches();

    match run(&matches) {
        Ok(()) => (),
        Err(e) => eprintln!("{}", e),
    }
}

fn run(args: &ArgMatches) -> Result<(), String> {
    let genes_path = args.value_of("genes").expect("should have genes arg");
    let annos_path = args.value_of("annotations").expect("should have annotations arg");
    let genes_out = args.value_of("genes_out").expect("should get genes output path");
    let annos_out = args.value_of("annotations_out").expect("should get annotations output path");
    let query = args.value_of("query").expect("should get query value");
    let segments = args.values_of("segment").expect("aspect should have some values");

    let segments: Vec<Segment> = segments.into_iter().map(|segment| {
        let split: Vec<&str> = segment.split(",").collect();
        let segment = (split[0], split[1]);
        Segment::try_from(segment).expect("should convert segment arg to Segment")
    }).collect();

    let mut genes_file = std::fs::File::open(genes_path)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut gene_reader = MetadataReader::new(BufReader::new(&mut genes_file));
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}",e ))?;
    let gene_metadata = gene_reader.metadata().expect("should capture gene metadata");
    let gene_headers = gene_reader.header().expect("should get gene headers");

    let mut annos_file = std::fs::File::open(annos_path)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata");
    let anno_headers = anno_reader.header().expect("should capture annotation header");

    let genes: Vec<Gene> = gene_records.iter()
        .map(|record| Gene::from_record(record))
        .collect();

    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_records.iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();

    let index: Index = Index::new(&genes, &annotations);
    let query = match query {
        "union" => Query::Union(segments),
        // "intersection" => Query::Intersection(segments),
        "intersection" => return Err(format!("Intersection queries are not yet implemented!")),
        _ => unreachable!(),
    };

    eprintln!("Executing query: {:?}", query);
    let result = query.execute(&index);

    let mut genes_out = std::fs::File::create(genes_out)
        .map_err(|e| format!("failed to create genes output file: {:?}", e))?;
    let mut genes_exporter = GafExporter::new(
        gene_metadata.to_string(),
        gene_headers.to_string(),
        result.genes_iter().map(|gene| gene.record));
    genes_exporter.write_all(&mut genes_out).expect("should write genes file");

    let mut annotations_out = std::fs::File::create(annos_out)
        .map_err(|e| format!("failed to create annotations output file: {:?}", e))?;
    let mut annotations_exporter = GafExporter::new(
        anno_metadata.to_string(),
        anno_headers.to_string(),
        result.annotations_iter().map(|anno| anno.record));
    annotations_exporter.write_all(&mut annotations_out)
        .map_err(|e| format!("failed to export data as GAF: {:?}", e))?;

    Ok(())
}
