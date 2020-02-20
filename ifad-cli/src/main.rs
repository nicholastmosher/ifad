use clap::{App, Arg, ArgMatches};
use std::io::{Write, BufReader};
use ifad::{MetadataReader, Annotation, Gene, Index, Aspect, AnnotationStatus, Segment};

fn app<'a, 'b>() -> clap::App<'a, 'b> {
    App::new("ifad")
        .arg(Arg::with_name("genes")
            .long("--genes")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations")
            .long("--annotations")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("output")
            .long("--output")
            .short("o")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("aspect")
            .long("--aspect")
            .possible_values(&["F", "P", "C"])
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("annotation_status")
            .long("--status")
            .possible_values(&["EXP", "OTHER", "UNKNOWN", "UNANNOTATED"])
            .required(true)
            .takes_value(true))
}

fn main() {
    let matches = app().get_matches();

    match run(&matches) {
        Ok(()) => (),
        Err(e) => eprintln!("{:?}", e),
    }
}

fn run(args: &ArgMatches) -> Result<(), String> {
    let genes_path = args.value_of("genes").expect("should have genes arg");
    let annos_path = args.value_of("annotations").expect("should have annotations arg");
    let out_path = args.value_of("output").expect("should have output arg");
    let aspect = args.value_of("aspect").expect("should have aspect");
    let status = args.value_of("annotation_status").expect("should have annotation status");

    let aspect = match aspect {
        "F" => Aspect::MolecularFunction,
        "P" => Aspect::BiologicalProcess,
        "C" => Aspect::CellularComponent,
        _ => unreachable!(),
    };

    let status = match status {
        "EXP" => AnnotationStatus::KnownExperimental,
        "OTHER" => AnnotationStatus::KnownOther,
        "UNKNOWN" => AnnotationStatus::Unknown,
        "UNANNOTATED" => AnnotationStatus::Unannotated,
        _ => unreachable!(),
    };

    let mut genes_file = std::fs::File::open(genes_path)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut gene_reader = MetadataReader::new(BufReader::new(&mut genes_file));
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}",e ))?;
    let _gene_metadata = gene_reader.metadata().expect("should capture gene metadata");

    let mut annos_file = std::fs::File::open(annos_path)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let _anno_metadata = anno_reader.metadata().expect("should capture annotation metadata");

    let genes: Vec<Gene> = gene_records.into_iter()
        .map(|record| Gene::from(record))
        .collect();

    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_records.into_iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();

    let index: Index = Index::new(&genes, &annotations);
    let segment = Segment::new(aspect, status);
    let result = segment.query(&index);

    let mut out_file = std::fs::File::create(out_path).expect("should open output file");
    write!(&mut out_file, "{:#?}", result)
        .map_err(|e| format!("failed to write to file: {:?}", e))?;

    Ok(())
}
