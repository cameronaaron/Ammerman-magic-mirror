"""Typer + Rich command-line interface.

    magic-mirror                     # live mirror (default), YuNet face tracking
    magic-mirror live --fullscreen
    magic-mirror poster out.png      # still image, no webcam needed
    magic-mirror list                # pretty table of your notable traits
    magic-mirror devices             # what removable-media genomes are visible
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table as RichTable

from .config import Settings, Source

app = typer.Typer(
    add_completion=False,
    no_args_is_help=False,
    help="Genetic Reflexions — a magic mirror that shows your genome, not your face.",
)
console = Console()
err = Console(stderr=True)


def _settings(ctx: typer.Context) -> Settings:
    return ctx.obj


def _log(message: str) -> None:
    err.print(message, style="dim")


# --------------------------------------------------------------------------- #
# Root: shared options, defaults to the live mirror
# --------------------------------------------------------------------------- #
@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    genome: Path | None = typer.Option(None, help="Path to a 23andMe raw file."),
    source: Source | None = typer.Option(None, help="Trait data source."),
    min_magnitude: float | None = typer.Option(
        None, "--min-magnitude", "-m", help="Minimum SNPedia magnitude to show."
    ),
    refresh: bool = typer.Option(False, help="Ignore the SNPedia cache; re-fetch."),
    snpedia: Path | None = typer.Option(None, help="CSV snapshot (fallback)."),
    device: int | None = typer.Option(None, help="Camera device id."),
    no_mirror: bool = typer.Option(False, "--no-mirror", help="Don't flip the image."),
    no_usb: bool = typer.Option(False, "--no-usb", help="Don't auto-detect USB genome."),
) -> None:
    overrides: dict = {}
    if genome is not None:
        overrides["genome_path"] = genome
    if source is not None:
        overrides["source"] = source
    if min_magnitude is not None:
        overrides["min_magnitude"] = min_magnitude
    if refresh:
        overrides["force_refresh"] = True
    if snpedia is not None:
        overrides["snpedia_path"] = snpedia
    if device is not None:
        overrides["device_id"] = device
    if no_mirror:
        overrides["mirror"] = False
    if no_usb:
        overrides["genome_from_usb"] = False
    ctx.obj = Settings(**overrides)

    if ctx.invoked_subcommand is None:
        _run_live(ctx.obj, fullscreen=False)


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #
@app.command()
def live(
    ctx: typer.Context,
    fullscreen: bool = typer.Option(False, help="Run the window fullscreen."),
) -> None:
    """Run the live webcam mirror."""
    _run_live(_settings(ctx), fullscreen=fullscreen)


def _run_live(settings: Settings, *, fullscreen: bool) -> None:
    from .app import run_live

    try:
        run_live(settings, fullscreen=fullscreen, log=_log)
    except (FileNotFoundError, ValueError, RuntimeError, ImportError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc


@app.command()
def poster(
    ctx: typer.Context,
    output: Path = typer.Argument(..., help="Where to write the PNG."),
    width: int = typer.Option(1280),
    height: int = typer.Option(800),
    at: float = typer.Option(1.2, help="Animation time in seconds."),
    background: Path | None = typer.Option(None, help="Optional background image."),
) -> None:
    """Render a still poster image (no webcam / OpenCV required)."""
    from .app import run_poster

    try:
        path = run_poster(
            _settings(ctx),
            output,
            size=(width, height),
            at_seconds=at,
            background=background,
            log=_log,
        )
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc
    console.print(f"Poster written → [green]{path}[/green]")


@app.command()
def share(
    ctx: typer.Context,
    photo: Path = typer.Argument(..., help="Existing photo to decorate."),
    output: Path = typer.Argument(..., help="Where to write the shareable PNG."),
    format: str = typer.Option(
        "portrait", "--format", "-f",
        help="Output size: square · portrait · story · wide.",
    ),
    no_detect: bool = typer.Option(
        False, "--no-detect", help="Skip face detection; centre the overlay."
    ),
    good: bool = typer.Option(
        False, "--good", help="Only show flattering 'good news' traits."
    ),
) -> None:
    """Turn an existing photo + your genome into an Instagram-ready image."""
    from .app import run_share

    try:
        path = run_share(
            _settings(ctx), photo, output,
            fmt=format, detect=not no_detect, positive_only=good, log=_log,
        )
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc
    console.print(f"Shareable image → [green]{path}[/green]  ([cyan]{format}[/cyan])")


@app.command()
def report(
    ctx: typer.Context,
    output: Path = typer.Argument("report.html", help="Where to write the HTML report."),
    csv_out: Path | None = typer.Option(
        None, "--csv", help="Also write the data as a flat CSV."
    ),
    no_genes: bool = typer.Option(
        False, "--no-genes", help="Skip fetching gene names from SNPedia."
    ),
) -> None:
    """Generate an interactive Promethease-style HTML report of every match."""
    from . import report as report_mod
    from .app import build_all_traits

    settings = _settings(ctx)
    try:
        genome, traits = build_all_traits(settings, log=_log)
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc

    genes: dict[str, str] = {}
    if not no_genes:
        from .snpedia.genes import fetch_genes

        genes = fetch_genes(
            (t.rsid for t in traits),
            cache_dir=settings.cache_dir,
            timeout=settings.http_timeout,
            log=_log,
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        report_mod.render_html(
            traits, genome=genome, genes=genes, title=settings.window_title
        ),
        encoding="utf-8",
    )
    counts = report_mod.report_counts(traits)
    console.print(
        f"Report → [green]{output}[/green]  "
        f"([bold]{counts['total']:,}[/bold] genotypes · "
        f"[green]{counts['Good']:,} good[/green] · "
        f"[red]{counts['Bad']:,} bad[/red] · "
        f"{counts['Neutral']:,} neutral)"
    )
    if csv_out is not None:
        path = report_mod.write_csv(traits, csv_out, genes)
        console.print(f"CSV → [green]{path}[/green]")


@app.command(name="gwas-report")
def gwas_report(
    ctx: typer.Context,
    output: Path = typer.Argument("gwas_report.html", help="Where to write the report."),
    csv_out: Path | None = typer.Option(None, "--csv", help="Also write a flat CSV."),
    limit: int = typer.Option(
        200, help="Max SNPs to query the GWAS Catalog for (bounds API calls)."
    ),
) -> None:
    """GWAS Catalog risk report: real odds ratios + risk-allele frequencies."""
    from . import report as report_mod
    from .app import build_all_traits
    from .snpedia.genes import fetch_genes
    from .snpedia.gwas import fetch_gwas

    settings = _settings(ctx)
    try:
        genome, all_traits = build_all_traits(settings, log=_log)
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc

    # The GWAS Catalog covers *common* variants; magnitude ~1.5–5 is the sweet
    # spot.  Very high magnitudes are rare Mendelian mutations it doesn't list,
    # so query the common range first when bounding the number of API calls.
    def _priority(t):
        return (0 if 1.5 <= t.magnitude < 5 else 1, -t.magnitude)

    traits = sorted(all_traits, key=_priority)[:limit]
    genotypes = {t.rsid: t.genotype for t in traits}
    rsids = list(genotypes)
    assocs = fetch_gwas(
        rsids, cache_dir=settings.cache_dir, timeout=settings.http_timeout, log=_log
    )
    genes = fetch_genes(
        rsids, cache_dir=settings.cache_dir, timeout=settings.http_timeout, log=_log
    )
    rows = report_mod.gwas_rows(genotypes, assocs, genes)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        report_mod.render_gwas_html(rows, genome=genome, title=settings.window_title),
        encoding="utf-8",
    )
    carried = sum(1 for r in rows if r["dosage"] > 0)
    console.print(
        f"GWAS report → [green]{output}[/green]  "
        f"([bold]{len(rows):,}[/bold] associations across {len(assocs):,} SNPs · "
        f"[yellow]{carried:,}[/yellow] risk alleles you carry)"
    )
    if csv_out is not None:
        path = report_mod.write_gwas_csv(rows, csv_out)
        console.print(f"CSV → [green]{path}[/green]")


@app.command(name="list")
def list_traits(ctx: typer.Context) -> None:
    """Print your notable traits as a table."""
    from .app import build_reflection

    try:
        reflection = build_reflection(_settings(ctx), log=_log)
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc

    table = RichTable(title=reflection.subheader, title_style="bold")
    table.add_column("mag", justify="right", style="bold")
    table.add_column("rsid", style="cyan")
    table.add_column("genotype", style="magenta")
    table.add_column("category", style="green")
    table.add_column("summary")
    for t in reflection.traits:
        table.add_row(
            f"{t.magnitude:g}", t.rsid, t.display_genotype, t.category, t.summary
        )
    console.print(table)


@app.command()
def chromosomes(
    ctx: typer.Context,
    compare: Path | None = typer.Option(
        None, help="Second genome to compare against (concordant/discordant SNPs)."
    ),
    output: Path | None = typer.Option(
        None, help="Write the table (or discordant SNPs, with --compare) to CSV."
    ),
    no_match: bool = typer.Option(
        False, "--no-match", help="Skip the SNPedia match column (faster, offline)."
    ),
) -> None:
    """Split the genome by chromosome; optionally diff it against another genome."""
    from . import chromosomes as chrom
    from .app import resolve_genome_path

    settings = _settings(ctx)
    try:
        genome_path = resolve_genome_path(settings, log=_log)
        if compare is not None:
            _print_comparison(chrom, genome_path, compare, output)
        else:
            _print_summary(chrom, settings, genome_path, output, no_match=no_match)
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc


def _print_summary(chrom, settings, genome_path, output, *, no_match: bool) -> None:
    from . import snpedia
    from .app import build_resolver

    table = None if no_match else snpedia.resolve_table(settings, log=_log)
    resolver = build_resolver(settings)
    stats = chrom.summarize(genome_path, table, resolver=resolver)

    rich = RichTable(title=f"Markers by chromosome — {genome_path.name}")
    rich.add_column("chr", style="cyan")
    rich.add_column("markers", justify="right")
    rich.add_column("diploid", justify="right")
    if not no_match:
        rich.add_column("matched", justify="right", style="green")
    for s in stats:
        row = [s.chromosome, f"{s.markers:,}", f"{s.diploid:,}"]
        if not no_match:
            row.append(f"{s.matched:,}")
        rich.add_row(*row)
    totals = [
        "Σ",
        f"{sum(s.markers for s in stats):,}",
        f"{sum(s.diploid for s in stats):,}",
    ]
    if not no_match:
        totals.append(f"{sum(s.matched for s in stats):,}")
    rich.add_section()
    rich.add_row(*totals, style="bold")
    console.print(rich)

    if output is not None:
        path = chrom.write_summary_csv(stats, output)
        console.print(f"Wrote summary → [green]{path}[/green]")


def _print_comparison(chrom, genome_path, compare, output) -> None:
    result = chrom.compare_genomes(genome_path, compare)
    rich = RichTable(
        title=f"{genome_path.name} vs {Path(compare).name} — shared SNPs by chromosome"
    )
    rich.add_column("chr", style="cyan")
    rich.add_column("shared", justify="right")
    rich.add_column("concordant", justify="right", style="green")
    rich.add_column("discordant", justify="right", style="yellow")
    for chromosome, (shared, conc, disc) in result.by_chromosome.items():
        rich.add_row(chromosome, f"{shared:,}", f"{conc:,}", f"{disc:,}")
    rich.add_section()
    rich.add_row(
        "Σ",
        f"{result.shared:,}",
        f"{result.concordant:,}",
        f"{result.discordant:,}",
        style="bold",
    )
    console.print(rich)
    console.print(f"Concordance: [bold]{result.concordance_rate:.2%}[/bold]")

    if output is not None:
        path = chrom.write_discordant_csv(result, output)
        console.print(
            f"Wrote {result.discordant:,} discordant SNPs → [green]{path}[/green]"
        )


@app.command()
def resolve(
    ctx: typer.Context,
    show: int = typer.Option(0, help="List up to N still-unresolved internal ids."),
) -> None:
    """Report how many 23andMe internal (i) ids translate to standard RSIDs."""
    from .app import build_resolver, resolve_genome_path
    from .genome import load_genome
    from .iids import IidResolver, is_internal_id

    settings = _settings(ctx)
    resolver = build_resolver(settings) or IidResolver()
    try:
        genome = load_genome(resolve_genome_path(settings, log=_log), resolver=resolver)
    except (FileNotFoundError, ValueError) as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc

    internal = [m for m in genome.calls if is_internal_id(m)]
    pct = (genome.resolved_ids / genome.internal_ids * 100) if genome.internal_ids else 0
    console.print(
        f"[bold]{genome.internal_ids:,}[/bold] internal ids · "
        f"[green]{genome.resolved_ids:,}[/green] resolved "
        f"([green]{pct:.1f}%[/green]) · "
        f"[yellow]{genome.internal_ids - genome.resolved_ids:,}[/yellow] unresolved"
    )
    console.print(f"Crosswalk: {settings.iid_map_path} ({len(resolver):,} mappings)")
    if show > 0:
        unresolved = resolver.unresolved(internal)
        for iid in unresolved[:show]:
            console.print(f"  [yellow]{iid}[/yellow]")
        if len(unresolved) > show:
            console.print(f"  … and {len(unresolved) - show:,} more")


@app.command()
def serve(
    ctx: typer.Context,
    host: str | None = typer.Option(None, help="Bind host (default 127.0.0.1)."),
    port: int | None = typer.Option(None, help="Bind port (default 8000)."),
) -> None:
    """Run the web upload/accounts server (needs the optional web extra)."""
    settings = _settings(ctx)
    try:
        import uvicorn

        from .web import create_app
    except (ImportError, RuntimeError) as exc:
        err.print(
            f"[red]error:[/red] {exc}\n"
            "Install the web extras: [cyan]pip install -e '.[web]'[/cyan]"
        )
        raise typer.Exit(1) from exc

    app_ = create_app(settings, data_dir=settings.web_data_dir)
    bind_host = host or settings.web_host
    bind_port = port or settings.web_port
    console.print(
        f"Serving Genetic Reflexions on [green]http://{bind_host}:{bind_port}[/green] "
        f"(data → {settings.web_data_dir})"
    )
    uvicorn.run(app_, host=bind_host, port=bind_port)


@app.command()
def snapshot(
    ctx: typer.Context,
    output: Path = typer.Argument(None, help="Where to write (default: result.csv)."),
) -> None:
    """Rebuild the offline SNPedia snapshot (result.csv) from a fresh live walk."""
    from . import snpedia
    from .snpedia import SnpediaUnavailable

    settings = _settings(ctx)
    target = output or settings.snpedia_path
    try:
        entries = snpedia.fetch_snapshot(timeout=settings.http_timeout, log=_log)
    except SnpediaUnavailable as exc:
        err.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(1) from exc

    path = snpedia.write_snapshot_csv(entries, target)
    good = sum(1 for e in entries if e.repute == "Good")
    bad = sum(1 for e in entries if e.repute == "Bad")
    console.print(
        f"Snapshot → [green]{path}[/green]  "
        f"([bold]{len(entries):,}[/bold] genotypes · "
        f"[green]{good:,} good[/green] · [red]{bad:,} bad[/red])"
    )


@app.command()
def devices(ctx: typer.Context) -> None:
    """List genome files found on removable media and in the local drop folder."""
    from .sources import iter_candidate_genomes, iter_folder_genomes

    settings = _settings(ctx)
    found = [("USB", p) for p in iter_candidate_genomes()]
    found += [("local", p) for p in iter_folder_genomes(settings.genome_inbox)]
    if not found:
        console.print(
            "No 23andMe genome files found on removable media or in "
            f"{settings.genome_inbox}."
        )
        return
    table = RichTable(title="Discovered genomes")
    table.add_column("where", style="green")
    table.add_column("path", style="cyan")
    for where, path in found:
        table.add_row(where, str(path))
    console.print(table)


if __name__ == "__main__":  # pragma: no cover
    app()
