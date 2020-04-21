use actix_web::{Scope, web};

pub mod genes;

pub fn routes(app: Scope) -> Scope {
    app
        .service(web::resource("genes")
            .route(web::get().to(genes::read)))
}
