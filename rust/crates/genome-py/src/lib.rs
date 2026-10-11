//! PyO3 bindings: expose `genome-core`'s fast parser to Python as `genome_core`.
//!
//! `magic_mirror.genome` imports this when present and falls back to the pure
//! Python parser otherwise, so the package works with or without the native
//! extension built.

use std::collections::HashMap;

use gcore::{build_traits, parse_genome, BuildOpts, SnpEntry, SnpTable};
use pyo3::prelude::*;

type ParsedTuple = (HashMap<String, String>, usize, usize, usize, usize);
// (rsid, genotype, magnitude, summary, category, repute)
type TraitTuple = (String, String, f64, String, String, String);
// (rsid, genotype, magnitude, repute, summary) as passed from Python's table.
type SnpRow = (String, String, f64, String, String);

/// Parse a 23andMe/AncestryDNA file at `path`, resolving internal ids via
/// `iid_map`.  Returns `(calls, total, diploid, internal, resolved)`.
#[pyfunction]
fn parse_genome_file(path: &str, iid_map: HashMap<String, String>) -> PyResult<ParsedTuple> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyOSError::new_err(e.to_string()))?;
    let g = parse_genome(&text, &iid_map);
    Ok((g.calls, g.total, g.diploid, g.internal, g.resolved))
}

/// Parse a genome **and** join it against a SNPedia table in one pass, returning
/// only the matched traits — the fast path (only a few hundred rows cross FFI,
/// not the whole 600k-marker dict).
#[pyfunction]
#[pyo3(signature = (path, iid_map, snp_rows, min_magnitude=0.0, require_summary=true, positive_only=false, limit=None))]
fn match_traits_file(
    path: &str,
    iid_map: HashMap<String, String>,
    snp_rows: Vec<SnpRow>,
    min_magnitude: f64,
    require_summary: bool,
    positive_only: bool,
    limit: Option<usize>,
) -> PyResult<Vec<TraitTuple>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyOSError::new_err(e.to_string()))?;
    let genome = parse_genome(&text, &iid_map);

    let mut table = SnpTable::with_capacity(snp_rows.len());
    for (rsid, genotype, magnitude, repute, summary) in snp_rows {
        table.insert((rsid, genotype), SnpEntry { magnitude, repute, summary });
    }
    let opts = BuildOpts { min_magnitude, require_summary, positive_only, limit };
    let traits = build_traits(&genome.calls, &table, &opts)
        .into_iter()
        .map(|t| (t.rsid, t.genotype, t.magnitude, t.summary, t.category, t.repute))
        .collect();
    Ok(traits)
}

#[pymodule]
fn genome_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_genome_file, m)?)?;
    m.add_function(wrap_pyfunction!(match_traits_file, m)?)?;
    m.add("__doc__", "Native genome parser + SNPedia join (Rust) for magic_mirror.")?;
    Ok(())
}
