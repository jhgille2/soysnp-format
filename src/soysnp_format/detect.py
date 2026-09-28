"""Detect which allele encoding a genotype file uses.

Detection works in two layers:

1. Syntactic evidence (needs no chip information). In A/B encoding the only
   allele letters that can ever appear are ``A`` and ``B``; the letter ``B``
   never occurs as a nucleotide. In FORWARD encoding genotypes are
   nucleotides, so ``C``, ``G`` or ``T`` can appear but ``B`` never does.
   Each SNP is therefore classified as AB-evidence (some genotype contains
   ``B``), FORWARD-evidence (some genotype contains ``C``/``G``/``T``),
   ambiguous (only ``A`` and/or missing calls observed), or conflicting
   (both kinds -- indicates a mixed-format file or data problems).

2. Lookup consistency (needs SNP ids matching a bundled chip table). For
   SNPs present in the lookup, genotypes whose alleles fit neither the A/B
   label set nor the SNP's forward allele pair are flagged as inconsistent
   with both encodings.

The verdict aggregates SNP-level evidence. A single SNP with a ``B`` call
proves A/B encoding for the file just as surely as a thousand do; the counts
are reported so the user can judge how much evidence the verdict rests on.
"""

from . import formats
from .formats import genotype_letters
from .lookup import load_lookup

AB = formats.AB
FORWARD = formats.FORWARD
INDETERMINATE = "INDETERMINATE"
MIXED = "MIXED"


class DetectionResult:
    """Outcome of :func:`detect_format`."""

    def __init__(self, verdict, confidence, evidence, declared_format=None,
                 chip_match=None, notes=()):
        self.verdict = verdict            # AB | FORWARD | MIXED | INDETERMINATE
        self.confidence = confidence      # high | moderate | low
        self.evidence = evidence          # dict of counts / lists
        self.declared_format = declared_format
        self.chip_match = chip_match      # {'chip': name, 'matched': n, 'total': n}
        self.notes = list(notes)

    def as_dict(self):
        return {
            "verdict": self.verdict,
            "confidence": self.confidence,
            "declared_format": self.declared_format,
            "chip_match": self.chip_match,
            "evidence": self.evidence,
            "notes": self.notes,
        }

    def summary(self):
        e = self.evidence
        lines = [
            f"Likely format: {self.verdict} (confidence: {self.confidence})",
            f"SNPs examined: {e['n_snps']}",
            f"  with A/B-only evidence (a 'B' call seen): {e['n_ab_evidence']}",
            f"  with forward-only evidence (C/G/T seen):  {e['n_forward_evidence']}",
            f"  ambiguous (only 'A' or no calls):         {e['n_ambiguous']}",
            f"  conflicting (both B and C/G/T seen):      {e['n_conflicting']}",
        ]
        if e.get("n_inconsistent"):
            lines.append(
                f"  inconsistent with both encodings:         {e['n_inconsistent']}"
            )
        if self.declared_format:
            lines.append(f"Format declared in file headers: {self.declared_format}")
            if self.declared_format in (AB, FORWARD) and self.verdict in (AB, FORWARD):
                agree = "agrees" if self.declared_format == self.verdict else "DISAGREES"
                lines.append(f"Header declaration {agree} with detected format.")
        if self.chip_match:
            m = self.chip_match
            lines.append(
                f"SNP ids matched against {m['chip']}: {m['matched']}/{m['total']}"
            )
        for note in self.notes:
            lines.append(f"Note: {note}")
        return "\n".join(lines)


def _classify_snp(letters, lookup_pair):
    """Classify one SNP's observed allele letters.

    Returns (evidence, inconsistent) where evidence is 'ab', 'forward',
    'ambiguous' or 'conflicting', and inconsistent flags genotypes fittable
    to neither encoding when a lookup pair is available.
    """
    has_b = "B" in letters
    has_cgt = bool(letters & {"C", "G", "T"})
    if has_b and has_cgt:
        evidence = "conflicting"
    elif has_b:
        evidence = "ab"
    elif has_cgt:
        evidence = "forward"
    else:
        evidence = "ambiguous"
    inconsistent = False
    if lookup_pair is not None and letters:
        a_fwd, b_fwd = lookup_pair
        fwd_pair = {a_fwd, b_fwd}
        # A genotype is explainable if its alleles are all A/B labels or all
        # forward nucleotides of this SNP.
        if not (letters <= {"A", "B"} or letters <= fwd_pair):
            inconsistent = True
    return evidence, inconsistent


def detect_format(data, lookup=None, chip="auto"):
    """Detect the allele encoding of parsed genotype ``data``.

    ``lookup`` is an optional :class:`~soysnp_format.lookup.ChipLookup` used
    for the consistency check and for chip matching. ``chip`` selects bundled
    tables when ``lookup`` is not given (``"auto"`` merges both chips).
    """
    if lookup is None:
        try:
            lookup = load_lookup(chip)
        except (FileNotFoundError, ValueError):
            lookup = None

    evidence_counts = {"ab": 0, "forward": 0, "ambiguous": 0, "conflicting": 0}
    inconsistent_snps = []
    matched = 0
    for snp in data.snps:
        letters = genotype_letters(gt for _, gt in data.iter_genotypes(snp))
        pair = lookup.alleles(snp) if lookup is not None and snp in lookup else None
        if pair is not None:
            matched += 1
        ev, inconsistent = _classify_snp(letters, pair)
        evidence_counts[ev] += 1
        if inconsistent:
            inconsistent_snps.append(snp)

    n = len(data.snps)
    n_ab = evidence_counts["ab"]
    n_fwd = evidence_counts["forward"]
    n_conf = evidence_counts["conflicting"]

    notes = []
    if n_conf or (n_ab and n_fwd):
        verdict = MIXED
        notes.append(
            f"{n_conf} SNP(s) show both 'B' calls and C/G/T calls within the SNP, "
            f"and {n_ab} SNP(s) look like A/B encoding while {n_fwd} look like "
            "forward encoding. The file mixes encodings, or SNP ids are "
            "misaligned with the data."
        )
    elif n_ab and not n_fwd:
        verdict = AB
    elif n_fwd and not n_ab:
        verdict = FORWARD
    else:
        verdict = INDETERMINATE
        notes.append(
            "No SNP shows a 'B' call (A/B evidence) or a C/G/T call (forward "
            "evidence); every genotype uses only the letter 'A' or is missing. "
            "Add more SNPs/samples, or specify --from explicitly when converting."
        )

    decisive = n_ab + n_fwd
    if verdict in (AB, FORWARD):
        if n_conf == 0 and decisive >= 10:
            confidence = "high"
        elif decisive >= 3:
            confidence = "moderate"
        else:
            confidence = "low"
            notes.append("Verdict rests on very few informative SNPs; treat with caution.")
    else:
        confidence = "low"

    declared = data.declared_format
    if declared and declared.startswith("UNKNOWN:"):
        notes.append(
            f"File headers name an encoding this tool does not translate "
            f"({declared.split(':', 1)[1]}); detection used the data itself."
        )
        declared = None
    if declared in (AB, FORWARD) and verdict in (AB, FORWARD) and declared != verdict:
        notes.append(
            f"File headers declare {declared} but the data look like {verdict}. "
            "The headers may be stale (e.g. edited by hand); the data were trusted."
        )

    chip_match = None
    if lookup is not None and n:
        chip_match = {"chip": lookup.name, "matched": matched, "total": n}
        if matched == 0:
            notes.append(
                "No SNP ids matched the bundled chip tables; detection used "
                "syntactic evidence only, and translation will need --lookup "
                "with a mapping for your chip."
            )

    evidence = {
        "n_snps": n,
        "n_ab_evidence": n_ab,
        "n_forward_evidence": n_fwd,
        "n_ambiguous": evidence_counts["ambiguous"],
        "n_conflicting": n_conf,
        "n_inconsistent": len(inconsistent_snps),
        "inconsistent_snps": inconsistent_snps[:25],
    }
    return DetectionResult(verdict, confidence, evidence,
                           declared_format=declared, chip_match=chip_match, notes=notes)
