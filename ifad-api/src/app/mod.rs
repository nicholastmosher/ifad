use std::sync::Arc;
use arc_swap::ArcSwap;
use actix_web::{App, HttpServer, web};
use ifad::Index;

pub mod v1;

pub struct Config {
    pub genes_file: String,
    pub annotations_file: String,
}

impl Config {
    pub fn from_env() -> Option<Config> {
        let genes_file = std::env::var("GENES_FILE").ok()?;
        let annotations_file = std::env::var("ANNOTATIONS_FILE").ok()?;
        Some(Config { genes_file, annotations_file })
    }
}

pub struct AppData {
    pub index: Arc<Index>,
    pub gene_metadata: String,
    pub gene_headers: String,
    pub anno_metadata: String,
    pub anno_headers: String,
}

pub async fn server(data: ArcSwap<AppData>) -> std::io::Result<()> {
    HttpServer::new(move || App::new()
        .data(data.clone())
        .configure(routes))
        .bind("127.0.0.1:8000")?
        .run()
        .await
}

fn routes(app: &mut web::ServiceConfig) {
    app.service(v1::routes(web::scope("/api/v1")));
}
