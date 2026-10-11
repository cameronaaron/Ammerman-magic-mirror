# Drop your 23andMe / AncestryDNA file here

The mirror looks for a genome in three places, in order:

1. an inserted **USB drive** (the installation's normal flow), then
2. **this folder** (`dna/`) — handy for testing or a fixed kiosk, then
3. the configured default (`sample/genome.txt`, a synthetic genome).

So to use your own data without a USB stick, just copy your raw-data export
(23andMe or AncestryDNA) into this folder:

```
cp ~/Downloads/AncestryDNA.txt dna/
magic-mirror list
```

Any file whose header looks like a 23andMe or AncestryDNA export (a vendor
banner, or an `rsid … genotype` / `rsid … allele1 allele2` column line) is
picked up automatically; the first match wins. Files you put here are
git-ignored, so your DNA never gets committed.
