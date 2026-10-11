//! Cloudflare Worker: genome-in, traits-out, at the edge.
//!
//! `POST /traits` with a 23andMe/AncestryDNA file as the request body returns
//! the SNPedia trait matches as JSON — the same `genome-core` parse + join the
//! Python app and the live mirror use, compiled to WASM.  The SNPedia snapshot
//! (`result.csv`) is bundled at build time, so the Worker is fully self-contained
//! and stateless (no D1/R2 needed for this read-only endpoint).
//!
//! Query params: `?min_magnitude=2&positive_only=true&limit=40`.

use std::collections::HashMap;
use std::sync::OnceLock;

use genome_core::{build_traits, parse_genome, parse_snp_csv, BuildOpts, SnpTable, Trait};
use serde::Serialize;
use worker::*;

// The public SNPedia snapshot, embedded at compile time.
static SNP_CSV: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../../result.csv"));

fn table() -> &'static SnpTable {
    static TABLE: OnceLock<SnpTable> = OnceLock::new();
    TABLE.get_or_init(|| parse_snp_csv(SNP_CSV))
}

#[derive(Serialize)]
struct TraitsResponse {
    total_markers: usize,
    diploid_markers: usize,
    matched: usize,
    traits: Vec<Trait>,
}

fn query_opts(url: &Url) -> BuildOpts {
    let mut opts = BuildOpts::default();
    for (k, v) in url.query_pairs() {
        match k.as_ref() {
            "min_magnitude" => opts.min_magnitude = v.parse().unwrap_or(0.0),
            "positive_only" => opts.positive_only = v == "true" || v == "1",
            "limit" => opts.limit = v.parse().ok(),
            _ => {}
        }
    }
    opts
}

#[event(fetch)]
async fn fetch(mut req: Request, _env: Env, _ctx: Context) -> Result<Response> {
    console_error_panic_hook::set_once();

    match (req.method(), req.path().as_str()) {
        (Method::Get, "/") => Response::ok(INDEX),
        (Method::Get, "/health") => Response::ok(format!("ok · {} genotypes", table().len())),
        (Method::Post, "/traits") => {
            let text = req.text().await?;
            if text.trim().is_empty() {
                return Response::error("empty body: POST a 23andMe/AncestryDNA file", 400);
            }
            let opts = query_opts(&req.url()?);
            let genome = parse_genome(&text, &HashMap::new());
            let traits = build_traits(&genome.calls, table(), &opts);
            let body = TraitsResponse {
                total_markers: genome.total,
                diploid_markers: genome.diploid,
                matched: traits.len(),
                traits,
            };
            let mut resp = Response::from_json(&body)?;
            resp.headers_mut().set("access-control-allow-origin", "*")?;
            Ok(resp)
        }
        _ => Response::error("Not Found", 404),
    }
}

const INDEX: &str = "Genetic Reflexions edge API\n\n\
POST /traits  (body: 23andMe/AncestryDNA raw file)\n\
  ?min_magnitude=2&positive_only=true&limit=40\n\
GET  /health\n";
