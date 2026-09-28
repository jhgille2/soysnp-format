/* soysnp.js — client-side port of the soysnp-format Python logic.
 *
 * Parses GenomeStudio exports and SNP matrices, detects the allele encoding
 * (Illumina A/B vs forward strand), and converts between encodings using a
 * per-SNP lookup table. No DOM dependencies: usable in node for testing.
 */

"use strict";

// ---------------------------------------------------------------------------
// formats
// ---------------------------------------------------------------------------

const MISSING_TOKENS = new Set([
  "", "-", "--", "---", "----", "?", ".", "./.",
  "NA", "N/A", "NAN", "NONE", "NULL",
  "NC", "NOCALL", "NO CALL", "NO_CALL",
  "NN", "N/N", "-/-", "--/--",
]);
const SEPARATORS = ["/", "|", " ", "\t", ":", ";", ","];
const AB = "AB";
const FORWARD = "FORWARD";
const INDETERMINATE = "INDETERMINATE";
const MIXED = "MIXED";

/** Parse one genotype token into [a1, a2], or null for missing calls. Throws on garbage. */
function normalizeGenotype(token) {
  if (token === null || token === undefined) return null;
  let text = String(token).trim().toUpperCase();
  if (MISSING_TOKENS.has(text)) return null;
  for (const sep of SEPARATORS) text = text.split(sep).join("");
  if (text.length !== 2) {
    throw new Error(`cannot parse genotype token ${JSON.stringify(token)}: expected two allele labels`);
  }
  const [a1, a2] = text;
  if (!"ACGTB".includes(a1) || !"ACGTB".includes(a2)) {
    throw new Error(`cannot parse genotype token ${JSON.stringify(token)}: unknown allele labels`);
  }
  return [a1, a2];
}

/** Distinct allele letters observed across genotypes (skips nulls). */
function genotypeLetters(genotypes) {
  const letters = new Set();
  for (const gt of genotypes) {
    if (gt === null || gt === undefined) continue;
    letters.add(gt[0]);
    letters.add(gt[1]);
  }
  return letters;
}

/** Render a genotype pair; alleles in alphabetical order for determinism. */
function formatGenotype(a1, a2, missing = "--") {
  if (a1 === null || a1 === undefined || a2 === null || a2 === undefined) return missing;
  return [a1, a2].sort().join("");
}

// ---------------------------------------------------------------------------
// delimited text helpers
// ---------------------------------------------------------------------------

function sniffDelimiter(sampleText) {
  const tabs = (sampleText.match(/\t/g) || []).length;
  const commas = (sampleText.match(/,/g) || []).length;
  const semis = (sampleText.match(/;/g) || []).length;
  if (tabs >= commas && tabs >= semis) return "\t";
  if (commas >= semis) return ",";
  return ";";
}

/** Minimal CSV parser: handles quoted fields and embedded newlines. */
function parseDelimited(text, delimiter) {
  const rows = [];
  let row = [], field = "", inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; }
        else inQuotes = false;
      } else field += c;
    } else if (c === '"') {
      inQuotes = true;
    } else if (c === delimiter) {
      row.push(field); field = "";
    } else if (c === "\r") {
      // skip; \n handles the break
    } else if (c === "\n") {
      row.push(field); field = "";
      rows.push(row); row = [];
    } else {
      field += c;
    }
  }
  if (field !== "" || row.length > 0) { row.push(field); rows.push(row); }
  return rows.filter(r => r.some(c => c.trim() !== ""));
}

// ---------------------------------------------------------------------------
// io: Final Report and wide matrices
// ---------------------------------------------------------------------------

const DECLARED_FORMAT_HINTS = {
  "AB": AB, "A/B": AB,
  "FORWARD": FORWARD, "PLUS": FORWARD,
  "TOP": null, "BOTTOM": null, "BOT": null, "MINUS": null, "DESIGN": null,
};

function declaredFormatFromHeader(columns) {
  for (const col of columns) {
    const m = /^\s*allele\s*[12]\s*-\s*(.+?)\s*$/i.exec(col);
    if (m) {
      const enc = m[1].trim().toUpperCase();
      return enc in DECLARED_FORMAT_HINTS ? DECLARED_FORMAT_HINTS[enc] : "UNKNOWN:" + enc;
    }
  }
  return null;
}

/** GenotypeData: { calls: Map(snp -> Map(sample -> [a1,a2]|null)), samples, snps, declaredFormat, layout, source } */
function makeGenotypeData(calls, samples, { declaredFormat = null, layout = null, source = "" } = {}) {
  return { calls, samples: [...samples], snps: [...calls.keys()], declaredFormat, layout, source };
}

/** Like normalizeGenotype, but returns null for garbage instead of throwing. */
function tryGenotype(token) {
  try { return normalizeGenotype(token); } catch { return null; }
}

/**
 * Parse a matrix-style Final Report [Data] section: the first data row holds
 * sample IDs (with an empty stub cell) and each following row holds one SNP
 * ID plus one genotype call per sample (e.g. AA/AB/BB). Returns a
 * GenotypeData, or null if the header row does not look like a sample-ID row.
 */
function readFinalReportMatrix(source, lines, dataStart, delim, columns) {
  if (columns.length < 2) return null;
  const stub = columns[0].trim().toLowerCase();
  const stubOk = !stub || ["snp", "snp name", "snpname", "marker", "locus", "sample", "sample id"].includes(stub);
  if (!stubOk) return null;
  const sampleCells = columns.slice(1);
  const samples = sampleCells.map(c => c.trim()).filter(Boolean);
  if (!samples.length) return null;
  const genoLike = sampleCells.filter(c => tryGenotype(c.trim()) !== null).length;
  if (genoLike > sampleCells.length / 2) return null;

  const dataRows = parseDelimited(lines.slice(dataStart + 1).join("\n"), delim);
  const calls = new Map();
  for (const row of dataRows) {
    if (!row.length || !row.some(c => c.trim())) continue;
    const snp = row[0].trim();
    if (!snp || snp.startsWith("[")) continue;
    if (!calls.has(snp)) calls.set(snp, new Map());
    row.slice(1).forEach((token, i) => {
      if (i < samples.length) calls.get(snp).set(samples[i], normalizeGenotype(token.trim()));
    });
  }
  if (!calls.size) return null;
  return makeGenotypeData(calls, samples, { declaredFormat: null, layout: "final-report", source });
}

function readFinalReport(text, source = "") {
  const lines = text.split(/\r?\n/);
  let dataStart = -1;
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim().toLowerCase() === "[data]") { dataStart = i + 1; break; }
  }
  if (dataStart < 0 || dataStart >= lines.length) {
    throw new Error("not a GenomeStudio Final Report ([Data] section not found)");
  }
  const delim = sniffDelimiter(lines.slice(dataStart, dataStart + 6).join("\n"));
  const columns = parseDelimited(lines[dataStart], delim)[0] || [];
  const colIndex = {};
  columns.forEach((c, i) => { colIndex[c.trim().toLowerCase()] = i; });
  const find = (...names) => { for (const n of names) if (n in colIndex) return colIndex[n]; return null; };

  const snpCol = find("snp name", "snpname", "snp", "locus");
  const sampleCol = find("sample id", "sampleid", "sample", "sample name", "samplename");
  let a1Col = find("allele1 - ab", "allele1-ab");
  let a2Col = find("allele2 - ab", "allele2-ab");
  if (a1Col === null) {
    a1Col = columns.findIndex(c => /^\s*allele\s*1\b/i.test(c));
    a2Col = columns.findIndex(c => /^\s*allele\s*2\b/i.test(c));
    if (a1Col < 0) a1Col = null;
    if (a2Col < 0) a2Col = null;
  }
  if (snpCol === null || sampleCol === null || a1Col === null || a2Col === null) {
    const matrix = readFinalReportMatrix(source, lines, dataStart, delim, columns);
    if (matrix) return matrix;
    throw new Error("Final Report header lacks SNP Name / Sample ID / Allele1 / Allele2 columns");
  }
  const declared = declaredFormatFromHeader(columns);

  const dataRows = parseDelimited(lines.slice(dataStart + 1).join("\n"), delim);
  const calls = new Map(), samples = [], seen = new Set();
  const need = Math.max(snpCol, sampleCol, a1Col, a2Col);
  for (const row of dataRows) {
    if (row.length <= need) continue;
    const snp = row[snpCol].trim(), sample = row[sampleCol].trim();
    if (!snp || !sample) continue;
    const gt = normalizeGenotype(row[a1Col].trim() + row[a2Col].trim());
    if (!seen.has(sample)) { seen.add(sample); samples.push(sample); }
    if (!calls.has(snp)) calls.set(snp, new Map());
    calls.get(snp).set(sample, gt);
  }
  return makeGenotypeData(calls, samples, { declaredFormat: declared, layout: "final-report", source });
}

function looksLikeSnpId(token, snpSet) {
  return snpSet.has(token.trim());
}

function readWideMatrix(text, snpSet = null, source = "") {
  const delim = sniffDelimiter(text.slice(0, 65536));
  let rows = parseDelimited(text, delim).map(r => r.map(c => c.trim()));
  // tolerate a GenomeStudio [Header]/[Data] preamble: parse from [Data]
  if (rows.length && rows[0].length && rows[0][0].startsWith("[")) {
    const di = rows.findIndex(r => r.length && r[0].trim().toLowerCase() === "[data]");
    if (di >= 0) rows = rows.slice(di + 1);
  }
  if (rows.length < 2 || rows[0].length < 2) {
    throw new Error("wide matrix needs at least 2 rows and 2 columns");
  }
  const header = rows[0];
  const firstCol = rows.slice(1).map(r => r[0]);
  const headerIds = header.slice(1).filter(Boolean);
  const firstColIds = firstCol.filter(Boolean);

  let snpsAreRows = null;
  if (snpSet && snpSet.size) {
    const headerHits = headerIds.filter(h => looksLikeSnpId(h, snpSet)).length;
    const colHits = firstColIds.filter(c => looksLikeSnpId(c, snpSet)).length;
    if (headerHits || colHits) snpsAreRows = colHits >= headerHits;
  }
  if (snpsAreRows === null) {
    const snpish = tok => /^(ss|rs|tm|sb)\d+$/i.test(tok);
    const headerHits = headerIds.filter(snpish).length;
    const colHits = firstColIds.filter(snpish).length;
    if (headerHits || colHits) snpsAreRows = colHits >= headerHits;
    else snpsAreRows = firstColIds.length >= headerIds.length;
  }

  const calls = new Map();
  let samples;
  if (snpsAreRows) {
    samples = headerIds;
    for (const row of rows.slice(1)) {
      const snp = row[0];
      if (!snp) continue;
      const m = new Map();
      headerIds.forEach((sample, i) => m.set(sample, normalizeGenotype(row[i + 1] ?? "")));
      calls.set(snp, m);
    }
  } else {
    samples = firstColIds;
    const snpIds = headerIds;
    for (const snp of snpIds) calls.set(snp, new Map());
    for (const row of rows.slice(1)) {
      const sample = row[0];
      if (!sample) continue;
      snpIds.forEach((snp, i) => calls.get(snp).set(sample, normalizeGenotype(row[i + 1] ?? "")));
    }
  }
  return makeGenotypeData(calls, samples, { layout: "wide", source });
}

function readGenotypes(text, layout = "auto", snpSet = null, source = "") {
  const head = text.slice(0, 4096).toLowerCase();
  const isFinalReport = head.includes("[data]") && head.includes("[header]");
  if (layout === "auto") layout = isFinalReport ? "final-report" : "wide";
  if (layout === "final-report") return readFinalReport(text, source);
  if (layout === "wide") return readWideMatrix(text, snpSet, source);
  throw new Error(`unknown layout ${JSON.stringify(layout)}`);
}

// ---------------------------------------------------------------------------
// detect
// ---------------------------------------------------------------------------

/** lookup: Map(snpId -> "AG") with both BARC and ss ids as keys, or null. */
function classifySnp(letters, lookupPair) {
  const hasB = letters.has("B");
  const hasCGT = [...letters].some(l => "CGT".includes(l));
  let evidence;
  if (hasB && hasCGT) evidence = "conflicting";
  else if (hasB) evidence = "ab";
  else if (hasCGT) evidence = "forward";
  else evidence = "ambiguous";
  let inconsistent = false;
  if (lookupPair && letters.size) {
    const fwdPair = new Set([lookupPair[0], lookupPair[1]]);
    const abOk = [...letters].every(l => l === "A" || l === "B");
    const fwdOk = [...letters].every(l => fwdPair.has(l));
    if (!abOk && !fwdOk) inconsistent = true;
  }
  return [evidence, inconsistent];
}

function detectFormat(data, lookup, chipName = "soysnp50k+soysnp6k") {
  const counts = { ab: 0, forward: 0, ambiguous: 0, conflicting: 0 };
  const inconsistentSnps = [];
  let matched = 0;
  for (const snp of data.snps) {
    const letters = genotypeLetters([...data.calls.get(snp).values()]);
    const pair = lookup && lookup.has(snp) ? lookup.get(snp) : null;
    if (pair) matched++;
    const [ev, inconsistent] = classifySnp(letters, pair);
    counts[ev]++;
    if (inconsistent) inconsistentSnps.push(snp);
  }
  const n = data.snps.length;
  const nAb = counts.ab, nFwd = counts.forward, nConf = counts.conflicting;
  const notes = [];
  let verdict;
  if (nConf || (nAb && nFwd)) {
    verdict = MIXED;
    notes.push(
      `${nConf} SNP(s) show both 'B' calls and C/G/T calls within the SNP, and ` +
      `${nAb} SNP(s) look like A/B encoding while ${nFwd} look like forward encoding. ` +
      "The file mixes encodings, or SNP ids are misaligned with the data."
    );
  } else if (nAb && !nFwd) verdict = AB;
  else if (nFwd && !nAb) verdict = FORWARD;
  else {
    verdict = INDETERMINATE;
    notes.push(
      "No SNP shows a 'B' call (A/B evidence) or a C/G/T call (forward evidence); " +
      "every genotype uses only the letter 'A' or is missing. Add more SNPs/samples, " +
      "or specify the source encoding explicitly when converting."
    );
  }
  const decisive = nAb + nFwd;
  let confidence;
  if (verdict === AB || verdict === FORWARD) {
    if (nConf === 0 && decisive >= 10) confidence = "high";
    else if (decisive >= 3) confidence = "moderate";
    else {
      confidence = "low";
      notes.push("Verdict rests on very few informative SNPs; treat with caution.");
    }
  } else confidence = "low";

  let declared = data.declaredFormat;
  if (declared && declared.startsWith("UNKNOWN:")) {
    notes.push(
      `File headers name an encoding this tool does not translate (${declared.slice(8)}); ` +
      "detection used the data itself."
    );
    declared = null;
  }
  if ((declared === AB || declared === FORWARD) &&
      (verdict === AB || verdict === FORWARD) && declared !== verdict) {
    notes.push(
      `File headers declare ${declared} but the data look like ${verdict}. ` +
      "The headers may be stale (e.g. edited by hand); the data were trusted."
    );
  }
  let chipMatch = null;
  if (lookup && n) {
    chipMatch = { chip: chipName, matched, total: n };
    if (matched === 0) {
      notes.push(
        "No SNP ids matched the bundled chip tables; detection used syntactic " +
        "evidence only, and translation needs a lookup table for your chip."
      );
    }
  }
  return {
    verdict, confidence, declaredFormat: declared, chipMatch, notes,
    evidence: {
      nSnps: n, nAbEvidence: nAb, nForwardEvidence: nFwd,
      nAmbiguous: counts.ambiguous, nConflicting: nConf,
      nInconsistent: inconsistentSnps.length,
      inconsistentSnps: inconsistentSnps.slice(0, 25),
    },
  };
}

// ---------------------------------------------------------------------------
// convert
// ---------------------------------------------------------------------------

class ConversionError extends Error {}

function alleleMap(fromFormat, toFormat, pair) {
  const [aFwd, bFwd] = [pair[0], pair[1]];
  if (fromFormat === AB && toFormat === FORWARD) return { A: aFwd, B: bFwd };
  if (fromFormat === FORWARD && toFormat === AB) {
    if (aFwd === bFwd) throw new ConversionError(`lookup row has identical forward alleles (${aFwd}); cannot invert`);
    const m = {};
    m[aFwd] = "A"; m[bFwd] = "B";
    return m;
  }
  throw new ConversionError(`unsupported conversion ${fromFormat} -> ${toFormat}`);
}

function convertGenotypes(data, fromFormat, toFormat, lookup,
                          { onUnknownSnp = "error", onBadAllele = "error" } = {}) {
  fromFormat = (fromFormat || "").toUpperCase();
  toFormat = (toFormat || "").toUpperCase();
  if (fromFormat === "AUTO") {
    const result = detectFormat(data, lookup);
    if (result.verdict === INDETERMINATE || result.verdict === MIXED) {
      throw new ConversionError(
        `cannot auto-detect source format (${result.verdict}, confidence ${result.confidence}); ` +
        "choose the source encoding explicitly."
      );
    }
    fromFormat = result.verdict;
  }
  if (fromFormat === toFormat) throw new ConversionError("source and target formats are identical; nothing to do");
  if (![AB, FORWARD].includes(fromFormat) || ![AB, FORWARD].includes(toFormat)) {
    throw new ConversionError(`unsupported conversion ${fromFormat} -> ${toFormat}; expected AB or FORWARD`);
  }
  if (!lookup) throw new ConversionError("no lookup table loaded");

  const unknownSnps = [];
  const badAlleles = [];
  const outCalls = new Map();
  for (const snp of data.snps) {
    if (!lookup.has(snp)) {
      unknownSnps.push(snp);
      if (onUnknownSnp === "skip") continue;
      continue;
    }
    const pair = lookup.get(snp);
    const amap = alleleMap(fromFormat, toFormat, pair);
    const allowed = fromFormat === AB ? new Set(["A", "B"]) : new Set([pair[0], pair[1]]);
    const m = new Map();
    for (const [sample, gt] of data.calls.get(snp)) {
      if (gt === null || gt === undefined) { m.set(sample, null); continue; }
      if (!allowed.has(gt[0]) || !allowed.has(gt[1])) {
        badAlleles.push([snp, sample, gt]);
        if (onBadAllele === "missing") { m.set(sample, null); continue; }
        continue;
      }
      m.set(sample, [amap[gt[0]], amap[gt[1]]]);
    }
    outCalls.set(snp, m);
  }
  const problems = [];
  if (unknownSnps.length && onUnknownSnp === "error") {
    problems.push(
      `${unknownSnps.length} SNP(s) not in the lookup table (e.g. ${unknownSnps.slice(0, 5).join(", ")}); ` +
      "these SNPs are not on the bundled chips"
    );
  }
  if (badAlleles.length && onBadAllele === "error") {
    const ex = badAlleles.slice(0, 5).map(([s, sm, gt]) => `${s}/${sm}=${gt.join("")}`);
    problems.push(
      `${badAlleles.length} genotype(s) use alleles outside the ${fromFormat} alphabet ` +
      `(e.g. ${ex.join(", ")})`
    );
  }
  if (problems.length) throw new ConversionError(problems.join("; "));
  return makeGenotypeData(outCalls, data.samples,
    { layout: "wide", source: `${data.source} [${fromFormat}->${toFormat}]` });
}

// ---------------------------------------------------------------------------
// writers
// ---------------------------------------------------------------------------

function writeWideMatrix(data, delimiter = "\t", missing = "--") {
  const lines = [["SNP", ...data.samples].join(delimiter)];
  for (const snp of data.snps) {
    const row = [snp];
    for (const sample of data.samples) {
      const gt = data.calls.get(snp).get(sample);
      row.push(formatGenotype(gt ? gt[0] : null, gt ? gt[1] : null, missing));
    }
    lines.push(row.join(delimiter));
  }
  return lines.join("\n") + "\n";
}

function writeLongTable(data, delimiter = "\t", missing = "--", encodingName = "Forward") {
  const lines = [[
    "SNP Name", "Sample ID", `Allele1 - ${encodingName}`, `Allele2 - ${encodingName}`,
  ].join(delimiter)];
  for (const snp of data.snps) {
    for (const sample of data.samples) {
      const gt = data.calls.get(snp).get(sample);
      const [a1, a2] = gt ? [...gt].sort() : [missing, missing];
      lines.push([snp, sample, a1, a2].join(delimiter));
    }
  }
  return lines.join("\n") + "\n";
}

// node / browser export
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    AB, FORWARD, INDETERMINATE, MIXED, ConversionError,
    normalizeGenotype, genotypeLetters, formatGenotype,
    sniffDelimiter, parseDelimited, declaredFormatFromHeader,
    readFinalReport, readWideMatrix, readGenotypes,
    classifySnp, detectFormat,
    convertGenotypes, writeWideMatrix, writeLongTable,
  };
}
