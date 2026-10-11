"""Generate a Promethease-style report from a genome × SNPedia match.

Produces a single **self-contained** HTML file (no CDN, no external fonts or
scripts — everything is inlined, so it opens offline and never phones home) that
replicates the core Promethease experience:

* every genotype SNPedia annotates, ranked by magnitude;
* **Repute** colour-coding — green "Good", red "Bad", slate "Neutral";
* live search over rsid / genotype / summary / category;
* repute toggles and a minimum-magnitude slider;
* good / bad / neutral tallies that update as you filter;
* deep links to each SNPedia page and genotype page.

Also emits the same data as a flat CSV for spreadsheet use.
"""

from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from .models import GenomeData, Trait
from .snpedia.gwas import GwasAssoc

_CSV_COLUMNS = ("rsid", "gene", "genotype", "magnitude", "repute", "category", "summary", "url")


def trait_to_row(trait: Trait, genes: dict[str, str] | None = None) -> dict:
    """A JSON/CSV-friendly view of one trait for the report."""
    return {
        "rsid": trait.rsid,
        "gene": (genes or {}).get(trait.rsid, ""),
        "genotype": trait.display_genotype,
        "magnitude": round(trait.magnitude, 3),
        "repute": trait.repute,
        "category": trait.category,
        "summary": trait.summary,
        "url": trait.snpedia_url,
        "gurl": trait.genotype_url,
    }


def report_counts(traits: list[Trait]) -> dict[str, int]:
    counts = {"Good": 0, "Bad": 0, "Neutral": 0}
    for trait in traits:
        counts[trait.repute] = counts.get(trait.repute, 0) + 1
    counts["total"] = len(traits)
    return counts


def write_csv(
    traits: list[Trait], path: str | Path, genes: dict[str, str] | None = None
) -> Path:
    """Write the report data as a flat CSV and return the path."""
    path = Path(path)
    genes = genes or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(_CSV_COLUMNS)
        for trait in traits:
            writer.writerow([
                trait.rsid, genes.get(trait.rsid, ""), trait.display_genotype,
                f"{trait.magnitude:g}", trait.repute, trait.category,
                trait.summary, trait.snpedia_url,
            ])
    return path


def render_html(
    traits: list[Trait],
    *,
    genome: GenomeData,
    genes: dict[str, str] | None = None,
    title: str = "Genetic Reflexions",
) -> str:
    """Return the full self-contained HTML report as a string."""
    ranked = sorted(traits, key=lambda t: (-t.magnitude, t.rsid))
    rows = [trait_to_row(t, genes) for t in ranked]
    counts = report_counts(ranked)

    data_json = json.dumps(rows, ensure_ascii=False).replace("</", "<\\/")
    subheader = (
        f"{genome.total_markers:,} markers read &nbsp;·&nbsp; "
        f"{len(ranked):,} SNPedia matches &nbsp;·&nbsp; "
        f"source: {html.escape(genome.source.name)}"
    )
    return (
        _TEMPLATE
        .replace("__TITLE__", html.escape(title))
        .replace("__SUBHEADER__", subheader)
        .replace("__GOOD__", str(counts["Good"]))
        .replace("__BAD__", str(counts["Bad"]))
        .replace("__NEUTRAL__", str(counts["Neutral"]))
        .replace("__TOTAL__", str(counts["total"]))
        .replace("__DATA__", data_json)
    )


# --------------------------------------------------------------------------- #
# GWAS risk report (odds ratios + risk-allele population frequencies)
# --------------------------------------------------------------------------- #
def gwas_rows(
    genotypes: dict[str, str],
    associations: dict[str, list[GwasAssoc]],
    genes: dict[str, str] | None = None,
) -> list[dict]:
    """Flatten GWAS associations into report rows with the user's risk-allele dose."""
    genes = genes or {}
    rows: list[dict] = []
    for rsid, assocs in associations.items():
        geno = genotypes.get(rsid, "")
        for a in assocs:
            rows.append({
                "rsid": rsid,
                "gene": genes.get(rsid, ""),
                "genotype": ";".join(geno) if geno else "—",
                "trait": a.trait,
                "risk_allele": a.risk_allele,
                "dosage": a.dosage(geno),
                "raf": round(a.risk_freq, 3) if a.risk_freq is not None else None,
                "pvalue": a.pvalue,
                "effect": a.effect,
                "url": f"https://www.snpedia.com/index.php/{rsid}",
            })
    # Most-carried, then most-significant first.
    rows.sort(key=lambda r: (-r["dosage"], r["pvalue"] if r["pvalue"] is not None else 1))
    return rows


def write_gwas_csv(rows: list[dict], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["rsid", "gene", "genotype", "trait", "risk_allele", "dosage",
            "raf", "pvalue", "effect"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(cols)
        for r in rows:
            writer.writerow([r[c] if r[c] is not None else "" for c in cols])
    return path


def render_gwas_html(
    rows: list[dict], *, genome: GenomeData, title: str = "Genetic Reflexions"
) -> str:
    carried = sum(1 for r in rows if r["dosage"] > 0)
    subheader = (
        f"{genome.total_markers:,} markers &nbsp;·&nbsp; "
        f"{len(rows):,} GWAS associations &nbsp;·&nbsp; "
        f"{carried:,} with a risk allele you carry &nbsp;·&nbsp; "
        f"source: GWAS Catalog (EBI)"
    )
    data_json = json.dumps(rows, ensure_ascii=False).replace("</", "<\\/")
    return (
        _GWAS_TEMPLATE
        .replace("__TITLE__", html.escape(title))
        .replace("__SUBHEADER__", subheader)
        .replace("__DATA__", data_json)
    )


# --------------------------------------------------------------------------- #
# Self-contained template.  Literal ``{``/``}`` are fine — we substitute with
# ``str.replace`` on __TOKENS__ rather than ``str.format``.
# --------------------------------------------------------------------------- #
_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ — genome report</title>
<style>
  :root{
    --bg:#101216; --panel:#171a21; --ink:#eef1f6; --muted:#9aa3b2;
    --line:rgba(255,255,255,.09); --good:#6ec6a8; --bad:#e2807a; --neutral:#9aa6ba;
    --goodbg:rgba(110,198,168,.10); --badbg:rgba(226,128,122,.10);
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
    font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;}
  a{color:#8ab4ff;text-decoration:none} a:hover{text-decoration:underline}
  header{position:sticky;top:0;z-index:10;background:rgba(16,18,22,.92);
    backdrop-filter:blur(8px);border-bottom:1px solid var(--line);padding:16px 22px;}
  h1{margin:0;font-size:20px;letter-spacing:.14em;text-transform:uppercase;font-weight:600}
  .sub{color:var(--muted);font-size:13px;margin-top:4px}
  .controls{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin-top:14px}
  #q{flex:1;min-width:220px;background:var(--panel);border:1px solid var(--line);
    color:var(--ink);border-radius:9px;padding:9px 12px;font:inherit}
  .chip{display:inline-flex;align-items:center;gap:7px;background:var(--panel);
    border:1px solid var(--line);border-radius:999px;padding:6px 12px;cursor:pointer;
    user-select:none;font-size:13px;color:var(--muted)}
  .chip input{accent-color:currentColor;margin:0}
  .chip.good{color:var(--good)} .chip.bad{color:var(--bad)} .chip.neutral{color:var(--neutral)}
  .chip.off{opacity:.4}
  .mag{display:flex;align-items:center;gap:8px;color:var(--muted);font-size:13px}
  .mag input{accent-color:#8ab4ff}
  main{max-width:1000px;margin:0 auto;padding:18px 22px 60px}
  .count{color:var(--muted);font-size:13px;margin:6px 0 14px}
  .card{border:1px solid var(--line);border-left:4px solid var(--neutral);
    background:var(--panel);border-radius:12px;padding:14px 16px;margin-bottom:10px}
  .card.Good{border-left-color:var(--good);background:linear-gradient(90deg,var(--goodbg),transparent 60%)}
  .card.Bad{border-left-color:var(--bad);background:linear-gradient(90deg,var(--badbg),transparent 60%)}
  .row1{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap}
  .glyph{font-weight:700;font-size:16px;width:1em;text-align:center}
  .Good .glyph{color:var(--good)} .Bad .glyph{color:var(--bad)} .Neutral .glyph{color:var(--neutral)}
  .rsid{font-family:ui-monospace,Menlo,monospace}
  .gene{font-size:11px;font-weight:700;letter-spacing:.04em;color:#c9b3ff;
    border:1px solid rgba(201,179,255,.35);border-radius:5px;padding:1px 6px}
  .geno{font-family:ui-monospace,Menlo,monospace;color:var(--muted)}
  .cat{margin-left:auto;font-size:11px;letter-spacing:.08em;text-transform:uppercase;
    color:var(--muted);border:1px solid var(--line);border-radius:999px;padding:2px 9px}
  .magbadge{font-weight:700;font-variant-numeric:tabular-nums;
    background:rgba(255,255,255,.06);border-radius:7px;padding:2px 9px;font-size:13px}
  .summary{margin:8px 0 0;color:var(--ink)}
  .empty{color:var(--muted);text-align:center;padding:50px 0}
  footer{color:var(--muted);font-size:12px;text-align:center;padding:24px;border-top:1px solid var(--line)}
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <div class="sub">__SUBHEADER__</div>
  <div class="controls">
    <input id="q" type="search" placeholder="Search rsid, genotype, condition, category…" autocomplete="off">
    <label class="chip good" data-rep="Good"><input type="checkbox" checked> Good <b id="cGood">__GOOD__</b></label>
    <label class="chip bad" data-rep="Bad"><input type="checkbox" checked> Bad <b id="cBad">__BAD__</b></label>
    <label class="chip neutral" data-rep="Neutral"><input type="checkbox" checked> Neutral <b id="cNeutral">__NEUTRAL__</b></label>
    <span class="mag">mag ≥ <input id="mag" type="range" min="0" max="10" step="0.5" value="0"><b id="magval">0</b></span>
  </div>
</header>
<main>
  <div class="count" id="count"></div>
  <div id="list"></div>
</main>
<footer>Data courtesy of <a href="https://www.snpedia.com/">SNPedia</a>.
For research and education only — not medical advice.</footer>
<script>
const DATA = __DATA__;
const list = document.getElementById('list');
const q = document.getElementById('q');
const mag = document.getElementById('mag');
const magval = document.getElementById('magval');
const repEnabled = {Good:true, Bad:true, Neutral:true};
const esc = s => s.replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const glyph = r => r==='Good' ? '+' : r==='Bad' ? '\\u2212' : '\\u2022';

function card(d){
  const summary = d.summary ? `<p class="summary">${esc(d.summary)}</p>` : '';
  const gene = d.gene ? `<span class="gene">${esc(d.gene)}</span>` : '';
  return `<div class="card ${d.repute}">
    <div class="row1">
      <span class="glyph">${glyph(d.repute)}</span>
      <a class="rsid" href="${d.url}" target="_blank" rel="noopener">${esc(d.rsid)}</a>
      ${gene}
      <a class="geno" href="${d.gurl}" target="_blank" rel="noopener">${esc(d.genotype)}</a>
      <span class="magbadge">${d.magnitude}</span>
      <span class="cat">${esc(d.category)}</span>
    </div>${summary}</div>`;
}

function render(){
  const terms = q.value.toLowerCase().split(/\\s+/).filter(Boolean);
  const floor = parseFloat(mag.value);
  const rows = DATA.filter(d => {
    if(!repEnabled[d.repute]) return false;
    if(d.magnitude < floor) return false;
    if(terms.length){
      const hay = (d.rsid+' '+d.gene+' '+d.genotype+' '+d.summary+' '+d.category).toLowerCase();
      if(!terms.some(t => hay.includes(t))) return false;
    }
    return true;
  });
  document.getElementById('count').textContent =
    rows.length.toLocaleString() + ' of ' + DATA.length.toLocaleString() + ' genotypes';
  list.innerHTML = rows.length
    ? rows.map(card).join('')
    : '<div class="empty">No matching genotypes.</div>';
}

q.addEventListener('input', render);
mag.addEventListener('input', () => { magval.textContent = mag.value; render(); });
document.querySelectorAll('.chip').forEach(chip => {
  const box = chip.querySelector('input');
  box.addEventListener('change', () => {
    repEnabled[chip.dataset.rep] = box.checked;
    chip.classList.toggle('off', !box.checked);
    render();
  });
});
render();
</script>
</body>
</html>
"""


_GWAS_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ — GWAS risk report</title>
<style>
  :root{--bg:#101216;--panel:#171a21;--ink:#eef1f6;--muted:#9aa3b2;
    --line:rgba(255,255,255,.09);--hot:#e2807a;--warm:#e0b877;--cool:#6ec6a8;}
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
    font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
  a{color:#8ab4ff;text-decoration:none} a:hover{text-decoration:underline}
  header{position:sticky;top:0;z-index:5;background:rgba(16,18,22,.94);
    backdrop-filter:blur(8px);border-bottom:1px solid var(--line);padding:16px 22px}
  h1{margin:0;font-size:20px;letter-spacing:.14em;text-transform:uppercase;font-weight:600}
  .sub{color:var(--muted);font-size:13px;margin-top:4px}
  #q{margin-top:12px;width:100%;max-width:420px;background:var(--panel);
    border:1px solid var(--line);color:var(--ink);border-radius:9px;padding:9px 12px;font:inherit}
  main{max-width:1100px;margin:0 auto;padding:14px 22px 60px;overflow-x:auto}
  table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
  th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
  th{color:var(--muted);font-weight:600;cursor:pointer;user-select:none;position:sticky;top:86px;background:var(--bg)}
  td.trait{white-space:normal;min-width:220px}
  .rsid{font-family:ui-monospace,Menlo,monospace}
  .gene{color:#c9b3ff;font-weight:700}
  .dose{font-weight:700;border-radius:6px;padding:1px 8px}
  .d2{background:rgba(226,128,122,.22);color:var(--hot)}
  .d1{background:rgba(224,184,119,.20);color:var(--warm)}
  .d0{background:rgba(110,198,168,.16);color:var(--cool)}
  footer{color:var(--muted);font-size:12px;text-align:center;padding:24px;border-top:1px solid var(--line)}
</style>
</head>
<body>
<header>
  <h1>__TITLE__ — GWAS</h1>
  <div class="sub">__SUBHEADER__</div>
  <input id="q" type="search" placeholder="Search trait, gene, rsid…" autocomplete="off">
</header>
<main><table id="t"><thead><tr>
  <th data-k="rsid">SNP</th><th data-k="gene">gene</th><th data-k="genotype">you</th>
  <th data-k="dosage">risk dose</th><th data-k="trait">GWAS trait</th>
  <th data-k="risk_allele">risk allele</th><th data-k="raf">pop. freq</th>
  <th data-k="pvalue">p-value</th>
</tr></thead><tbody id="body"></tbody></table></main>
<footer>Associations from the <a href="https://www.ebi.ac.uk/gwas/">GWAS Catalog</a> ·
frequencies are the reported risk-allele frequency. Research/education only.</footer>
<script>
const DATA = __DATA__;
const body = document.getElementById('body');
const q = document.getElementById('q');
let sortKey = 'dosage', asc = false;
const esc = s => (''+s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const fmtP = p => p==null ? '' : (p < 1e-3 ? p.toExponential(0) : (''+p));

function row(d){
  const dl = d.dosage===2?'d2':d.dosage===1?'d1':'d0';
  return `<tr>
    <td><a class="rsid" href="${d.url}" target="_blank" rel="noopener">${esc(d.rsid)}</a></td>
    <td class="gene">${esc(d.gene||'')}</td>
    <td class="rsid">${esc(d.genotype)}</td>
    <td><span class="dose ${dl}">${d.dosage}×${esc(d.risk_allele)}</span></td>
    <td class="trait">${esc(d.trait)}</td>
    <td class="rsid">${esc(d.risk_allele)}</td>
    <td>${d.raf==null?'':esc(d.raf)}</td>
    <td>${fmtP(d.pvalue)}</td></tr>`;
}
function render(){
  const terms = q.value.toLowerCase().split(/\\s+/).filter(Boolean);
  let rows = DATA.filter(d => !terms.length ||
    terms.some(t => (d.rsid+' '+d.gene+' '+d.trait).toLowerCase().includes(t)));
  rows.sort((a,b)=>{
    let x=a[sortKey], y=b[sortKey];
    if(x==null)x=asc?Infinity:-Infinity; if(y==null)y=asc?Infinity:-Infinity;
    if(typeof x==='string')return asc?(''+x).localeCompare(y):(''+y).localeCompare(x);
    return asc?x-y:y-x;
  });
  body.innerHTML = rows.map(row).join('') || '<tr><td colspan=8>No matches.</td></tr>';
}
q.addEventListener('input', render);
document.querySelectorAll('th').forEach(th=>th.addEventListener('click',()=>{
  const k=th.dataset.k; asc = k===sortKey ? !asc : false; sortKey=k; render();
}));
render();
</script>
</body>
</html>
"""
