/* app.js — UI wiring for the soysnp-format web tool. Requires soysnp.js. */
"use strict";

const $ = (id) => document.getElementById(id);

const EXAMPLE =
  "SNP\tsample1\tsample2\tsample3\n" +
  "BARC_1.01_Gm01_2033_G_A\tAA\tAB\tBB\n" +
  "BARC_1.01_Gm01_25990_C_T\tAB\tBB\t--\n" +
  "BARC_1.01_Gm01_29658_A_G\tAA\tAA\tAB\n";

let lookup = null;      // Map(snpId -> "AG")
let lookupLoading = null;
let data = null;        // parsed GenotypeData
let fileName = "";

async function ensureLookup() {
  if (lookup) return lookup;
  if (!lookupLoading) {
    lookupLoading = fetch("assets/lookups.json").then(async (resp) => {
      if (!resp.ok) throw new Error(`could not load lookup table (HTTP ${resp.status})`);
      const rows = await resp.json();
      const map = new Map();
      for (const [shortId, ssId, pair] of rows) {
        map.set("BARC_1.01_" + shortId, pair);
        if (ssId) map.set(ssId, pair);
      }
      lookup = map;
      return map;
    });
  }
  return lookupLoading;
}

function setStatus(msg) { $("fileStatus").textContent = msg; }

async function loadText(text, name) {
  fileName = name;
  try {
    const lk = await ensureLookup();
    data = readGenotypes(text, $("layout").value, lk, name);
    setStatus(`Loaded ${name}: ${data.snps.length} SNPs x ${data.samples.length} samples ` +
              `(${data.layout === "final-report" ? "GenomeStudio Final Report" : "wide matrix"}).`);
    $("detectBtn").disabled = false;
    $("convertBtn").disabled = false;
    $("detectOut").classList.add("hidden");
    $("preview").classList.add("hidden");
    $("download").classList.add("hidden");
    $("convertStatus").textContent = "";
  } catch (e) {
    data = null;
    $("detectBtn").disabled = true;
    $("convertBtn").disabled = true;
    setStatus(`Could not parse ${name}: ${e.message}`);
  }
}

// --- file input ---
const drop = $("drop"), fileInput = $("file");
drop.addEventListener("click", () => fileInput.click());
drop.addEventListener("dragover", (e) => { e.preventDefault(); drop.classList.add("over"); });
drop.addEventListener("dragleave", () => drop.classList.remove("over"));
drop.addEventListener("drop", (e) => {
  e.preventDefault(); drop.classList.remove("over");
  const f = e.dataTransfer.files[0];
  if (f) f.text().then((t) => loadText(t, f.name));
});
fileInput.addEventListener("change", () => {
  const f = fileInput.files[0];
  if (f) f.text().then((t) => loadText(t, f.name));
});
$("example").addEventListener("click", () => loadText(EXAMPLE, "example_ab.tsv"));

// --- detection ---
$("detectBtn").addEventListener("click", async () => {
  if (!data) return;
  const lk = await ensureLookup();
  let result;
  try {
    result = detectFormat(data, lk);
  } catch (e) {
    $("detectOut").classList.remove("hidden");
    $("verdict").textContent = "error";
    $("notes").innerHTML = `<div class="note">${escapeHtml(e.message)}</div>`;
    return;
  }
  $("detectOut").classList.remove("hidden");
  $("verdict").textContent = result.verdict;
  const conf = $("confidence");
  conf.textContent = `(confidence: ${result.confidence})`;
  conf.className = "conf-" + result.confidence;
  const ev = result.evidence;
  $("cN").textContent = ev.nSnps;
  $("cAb").textContent = ev.nAbEvidence;
  $("cFwd").textContent = ev.nForwardEvidence;
  $("cAmb").textContent = ev.nAmbiguous;
  $("cConf").textContent = ev.nConflicting;
  $("cInc").textContent = ev.nInconsistent +
    (ev.inconsistentSnps.length ? ` (${ev.inconsistentSnps.slice(0, 5).join(", ")}${ev.nInconsistent > 5 ? ", ..." : ""})` : "");
  $("declaredLine").textContent = result.declaredFormat
    ? `Format declared in file headers: ${result.declaredFormat}` : "";
  $("chipLine").textContent = result.chipMatch
    ? `SNP ids matched against ${result.chipMatch.chip}: ${result.chipMatch.matched}/${result.chipMatch.total}` : "";
  $("notes").innerHTML = result.notes.map((n) => `<div class="note">${escapeHtml(n)}</div>`).join("");
});

// --- conversion ---
$("convertBtn").addEventListener("click", async () => {
  if (!data) return;
  const lk = await ensureLookup();
  const from = $("fromFmt").value, to = $("toFmt").value;
  $("convertStatus").textContent = "Converting...";
  $("preview").classList.add("hidden");
  $("download").classList.add("hidden");
  try {
    const out = convertGenotypes(data, from, to, lk, {
      onUnknownSnp: $("skipUnknown").checked ? "skip" : "error",
      onBadAllele: $("badMissing").checked ? "missing" : "error",
    });
    const text = $("outLayout").value === "long"
      ? writeLongTable(out, "\t", "--", to === "AB" ? "AB" : "Forward")
      : writeWideMatrix(out);
    const blob = new Blob([text], { type: "text/tab-separated-values" });
    const url = URL.createObjectURL(blob);
    const a = $("download");
    a.href = url;
    a.download = fileName.replace(/\.[^.]*$/, "") + `_${to.toLowerCase()}.tsv`;
    a.classList.remove("hidden");
    $("preview").textContent = text.split("\n").slice(0, 12).join("\n") +
      (text.split("\n").length > 12 ? "\n..." : "");
    $("preview").classList.remove("hidden");
    $("convertStatus").textContent =
      `Converted ${out.snps.length} SNPs x ${out.samples.length} samples (${from.toUpperCase()} -> ${to}).`;
  } catch (e) {
    $("convertStatus").textContent = `Conversion failed: ${e.message}`;
  }
});

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}
