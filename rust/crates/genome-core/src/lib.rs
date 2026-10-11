//! Portable genome parsing + SNPedia trait-building.
//!
//! This is the single source of truth shared by the Python bindings
//! (`genome-py`, via PyO3) and the Cloudflare Worker (`genome-worker`, via
//! WASM).  It mirrors the pure-Python implementation in `magic_mirror/` so the
//! three run-times stay in lock-step: 23andMe **and** AncestryDNA parsing,
//! order-independent genotype matching, category inference, and SNPedia Repute
//! (Good / Bad / Neutral) with a category fallback.

use std::collections::HashMap;

use serde::Serialize;

// --------------------------------------------------------------------------- //
// Genotype helpers
// --------------------------------------------------------------------------- //
#[inline]
fn is_base(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

/// Strand/order-independent key for a diploid call, else `None` (`"AG"`→`"AG"`,
/// `"GA"`→`"AG"`; no-calls, indels and single alleles are rejected).
pub fn normalize_genotype(genotype: &str) -> Option<String> {
    let g = genotype.trim().to_ascii_uppercase();
    let b: Vec<char> = g.chars().collect();
    if b.len() == 2 && is_base(b[0]) && is_base(b[1]) {
        let (x, y) = if b[0] <= b[1] { (b[0], b[1]) } else { (b[1], b[0]) };
        return Some(format!("{x}{y}"));
    }
    None
}

/// A 23andMe internal id: `i` followed by digits.
pub fn is_internal_id(marker: &str) -> bool {
    let b = marker.as_bytes();
    b.len() >= 2 && (b[0] == b'i' || b[0] == b'I') && b[1..].iter().all(u8::is_ascii_digit)
}

fn is_rsid(marker: &str) -> bool {
    let b = marker.as_bytes();
    b.len() >= 3
        && (b[0] == b'r' || b[0] == b'R')
        && (b[1] == b's' || b[1] == b'S')
        && b[2..].iter().all(u8::is_ascii_digit)
}

// AncestryDNA encodes the sex chromosomes / mitochondria numerically.
fn map_ancestry_chrom(c: &str) -> &str {
    match c {
        "23" => "X",
        "24" => "Y",
        "25" => "XY",
        "26" => "MT",
        other => other,
    }
}

// --------------------------------------------------------------------------- //
// Genome parsing (23andMe + AncestryDNA)
// --------------------------------------------------------------------------- //
#[derive(Debug, Default)]
pub struct ParsedGenome {
    pub calls: HashMap<String, String>,
    pub total: usize,
    pub diploid: usize,
    pub internal: usize,
    pub resolved: usize,
}

/// Parse raw genome text into `rsid → normalized genotype` calls.
///
/// `iid_map` translates internal `i` ids to standard RSIDs (the original key is
/// kept as an alias too).  Handles both the 23andMe 4-column and AncestryDNA
/// 5-column layouts, told apart per row by column count.
pub fn parse_genome(text: &str, iid_map: &HashMap<String, String>) -> ParsedGenome {
    let mut g = ParsedGenome::default();
    for raw in text.lines() {
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').map(str::trim).collect();
        if fields[0].is_empty() || fields[0].eq_ignore_ascii_case("rsid") {
            continue;
        }
        let genotype = match fields.len() {
            4 => fields[3].to_string(),                        // 23andMe
            5 => format!("{}{}", fields[3], fields[4]),        // AncestryDNA
            _ => continue,
        };
        let _chrom = if fields.len() == 5 {
            map_ancestry_chrom(fields[1])
        } else {
            fields[1]
        };
        g.total += 1;
        let key = match normalize_genotype(&genotype) {
            Some(k) => k,
            None => continue,
        };
        g.diploid += 1;
        let marker = fields[0].to_ascii_lowercase();
        g.calls.entry(marker.clone()).or_insert_with(|| key.clone());
        if is_internal_id(&marker) {
            g.internal += 1;
            if let Some(canon) = iid_map.get(&marker) {
                if canon != &marker {
                    g.resolved += 1;
                    g.calls.entry(canon.clone()).or_insert(key);
                }
            }
        }
    }
    g
}

// --------------------------------------------------------------------------- //
// Category + repute (ports magic_mirror/traits.py)
// --------------------------------------------------------------------------- //
const CATEGORY_KEYWORDS: &[(&str, &[&str])] = &[
    ("longevity", &["longevity", "lifespan", "aging", "ageing", "nonagenarian"]),
    ("cognition", &["memory", "optimis", "empath", "stress", "intelligence", "iq"]),
    ("performance", &["muscle", "sprint", "power", "athlet", "endurance", "sleeper"]),
    ("pharma", &["drug", "plavix", "warfarin", "detoxif", "response to treatment"]),
    ("protection", &["reduced", "resistance", "decreased", "lower risk", "protect", "half the"]),
    ("risk", &["risk", "cancer", "disease", "asthma", "gout", "infarction", "sclerosis", "diabetes", "carrier"]),
    ("appearance", &["earwax", "hair", "eye", "skin", "odour", "odor", "height", "bald", "light-skinned", "myopia"]),
    ("metabolism", &["caffeine", "alcohol", "lactose", "metaboli", "diet", "weight", "cholesterol"]),
];

const POSITIVE_CATEGORIES: &[&str] =
    &["performance", "longevity", "cognition", "protection", "appearance"];

pub fn categorize(summary: &str) -> &'static str {
    let text = summary.to_ascii_lowercase();
    for (category, keywords) in CATEGORY_KEYWORDS {
        if keywords.iter().any(|k| text.contains(k)) {
            return category;
        }
    }
    "trait"
}

pub fn resolve_repute(raw: &str, category: &str) -> String {
    match raw.trim().to_ascii_lowercase().as_str() {
        "good" => "Good".to_string(),
        "bad" => "Bad".to_string(),
        _ => match category {
            "protection" => "Good".to_string(),
            "risk" => "Bad".to_string(),
            _ => "Neutral".to_string(),
        },
    }
}

pub fn is_positive(repute: &str, category: &str) -> bool {
    if repute == "Bad" {
        return false;
    }
    repute == "Good" || POSITIVE_CATEGORIES.contains(&category)
}

// --------------------------------------------------------------------------- //
// SNPedia table + trait building
// --------------------------------------------------------------------------- //
#[derive(Debug, Clone)]
pub struct SnpEntry {
    pub magnitude: f64,
    pub repute: String,
    pub summary: String,
}

pub type SnpTable = HashMap<(String, String), SnpEntry>;

fn parse_genotype_title(title: &str) -> Option<(String, String)> {
    let t = title.trim();
    let open = t.find('(')?;
    let id = t[..open].trim().to_ascii_lowercase();
    if !(is_rsid(&id) || is_internal_id(&id)) {
        return None;
    }
    let inner = t[open..].trim_start_matches('(').trim_end_matches(')');
    let alleles: Vec<&str> = inner.split(';').collect();
    if alleles.len() != 2 {
        return None;
    }
    let g = normalize_genotype(&format!("{}{}", alleles[0], alleles[1]))?;
    Some((id, g))
}

/// Parse a `result.csv`-style snapshot (`,Magnitude[,Repute],Summary`) into a table.
pub fn parse_snp_csv(text: &str) -> SnpTable {
    let mut table = SnpTable::new();
    let mut reader = csv::ReaderBuilder::new()
        .flexible(true)
        .from_reader(text.as_bytes());
    let headers = reader.headers().cloned().unwrap_or_default();
    let repute_idx = headers
        .iter()
        .position(|h| h.trim().eq_ignore_ascii_case("repute"));
    for record in reader.records().flatten() {
        let title = record.get(0).unwrap_or("");
        let (rsid, genotype) = match parse_genotype_title(title) {
            Some(v) => v,
            None => continue,
        };
        let magnitude = record.get(1).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
        let repute = repute_idx
            .and_then(|i| record.get(i))
            .unwrap_or("")
            .trim()
            .to_string();
        // Summary is the last column (index 2 without repute, 3 with).
        let summary_idx = if repute_idx.is_some() { 3 } else { 2 };
        let summary = record.get(summary_idx).unwrap_or("").trim().to_string();
        table.insert((rsid, genotype), SnpEntry { magnitude, repute, summary });
    }
    table
}

#[derive(Debug, Serialize)]
pub struct Trait {
    pub rsid: String,
    pub genotype: String,
    pub magnitude: f64,
    pub summary: String,
    pub category: String,
    pub repute: String,
}

#[derive(Debug, Clone)]
pub struct BuildOpts {
    pub min_magnitude: f64,
    pub require_summary: bool,
    pub positive_only: bool,
    pub limit: Option<usize>,
}

impl Default for BuildOpts {
    fn default() -> Self {
        Self { min_magnitude: 0.0, require_summary: true, positive_only: false, limit: None }
    }
}

/// Join genome calls against a SNPedia table, ranked most-notable first.
pub fn build_traits(
    calls: &HashMap<String, String>,
    table: &SnpTable,
    opts: &BuildOpts,
) -> Vec<Trait> {
    let mut traits = Vec::new();
    for (rsid, geno) in calls {
        let entry = match table.get(&(rsid.clone(), geno.clone())) {
            Some(e) => e,
            None => continue,
        };
        if entry.magnitude < opts.min_magnitude {
            continue;
        }
        if opts.require_summary && entry.summary.is_empty() {
            continue;
        }
        let category = categorize(&entry.summary);
        let repute = resolve_repute(&entry.repute, category);
        if opts.positive_only && !is_positive(&repute, category) {
            continue;
        }
        traits.push(Trait {
            rsid: rsid.clone(),
            genotype: format!("{};{}", &geno[0..1], &geno[1..2]),
            magnitude: entry.magnitude,
            summary: entry.summary.clone(),
            category: category.to_string(),
            repute,
        });
    }
    traits.sort_by(|a, b| {
        b.magnitude
            .partial_cmp(&a.magnitude)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.rsid.cmp(&b.rsid))
    });
    if let Some(limit) = opts.limit {
        traits.truncate(limit);
    }
    traits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_order_independent() {
        assert_eq!(normalize_genotype("GA").as_deref(), Some("AG"));
        assert_eq!(normalize_genotype("--"), None);
        assert_eq!(normalize_genotype("A"), None);
    }

    #[test]
    fn parses_both_vendor_formats() {
        let m = HashMap::new();
        let t23 = " rsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAG\ni9\t1\t2\tCC\n";
        let g = parse_genome(t23, &m);
        assert_eq!(g.diploid, 2);
        assert_eq!(g.calls["rs1"], "AG");
        assert_eq!(g.internal, 1);

        let anc = "rsid\tchromosome\tposition\tallele1\tallele2\nrs1\t1\t1\tG\tA\nrs2\t23\t2\t0\t0\n";
        let g = parse_genome(anc, &m);
        assert_eq!(g.calls["rs1"], "AG"); // split alleles joined + sorted
        assert_eq!(g.diploid, 1); // 0/0 no-call dropped
    }

    #[test]
    fn resolves_internal_ids() {
        let mut m = HashMap::new();
        m.insert("i9".to_string(), "rs900".to_string());
        let g = parse_genome("rs1\t1\t1\tAG\ni9\t1\t2\tCC\n", &m);
        assert_eq!(g.resolved, 1);
        assert_eq!(g.calls["rs900"], "CC");
    }

    #[test]
    fn builds_traits_with_repute() {
        let mut calls = HashMap::new();
        calls.insert("rs1".to_string(), "AG".to_string());
        let mut table = SnpTable::new();
        table.insert(
            ("rs1".to_string(), "AG".to_string()),
            SnpEntry { magnitude: 3.0, repute: "Good".into(), summary: "nice".into() },
        );
        let traits = build_traits(&calls, &table, &BuildOpts::default());
        assert_eq!(traits.len(), 1);
        assert_eq!(traits[0].genotype, "A;G");
        assert_eq!(traits[0].repute, "Good");
    }

    #[test]
    fn parses_snp_csv_with_repute() {
        let csv = ",Magnitude,Repute,Summary\nRs17822931(T;T),2.5,Good,\"Dry earwax, no odour\"\n";
        let table = parse_snp_csv(csv);
        let e = table.get(&("rs17822931".to_string(), "TT".to_string())).unwrap();
        assert_eq!(e.repute, "Good");
        assert_eq!(e.magnitude, 2.5);
    }
}
