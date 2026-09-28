import json
import subprocess
import sys

from soysnp_format.cli import main

from conftest import AB_WIDE, FORWARD_WIDE


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return str(p)


def _lookup_csv(tmp_path, lookup):
    p = tmp_path / "lookup.csv"
    lines = ["snp,a_forward,b_forward"]
    lines += [f"{snp},{a},{b}" for snp, (a, b) in sorted(lookup.mapping.items())]
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_cli_detect(tmp_path, lookup, capsys):
    inp = _write(tmp_path, "in.tsv", AB_WIDE)
    lkp = _lookup_csv(tmp_path, lookup)
    assert main(["detect", inp, "--lookup", lkp]) == 0
    out = capsys.readouterr().out
    assert "Likely format: AB" in out


def test_cli_detect_json(tmp_path, lookup, capsys):
    inp = _write(tmp_path, "in.tsv", FORWARD_WIDE)
    lkp = _lookup_csv(tmp_path, lookup)
    assert main(["detect", inp, "--lookup", lkp, "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["verdict"] == "FORWARD"


def test_cli_convert(tmp_path, lookup):
    inp = _write(tmp_path, "in.tsv", AB_WIDE)
    lkp = _lookup_csv(tmp_path, lookup)
    outp = str(tmp_path / "out.tsv")
    assert main(["convert", inp, "--lookup", lkp, "--from", "AB",
                 "--to", "FORWARD", "-o", outp]) == 0
    lines = open(outp).read().splitlines()
    assert lines[1] == "snp1\tAA\tAG\tGG"


def test_cli_convert_auto_from(tmp_path, lookup):
    inp = _write(tmp_path, "in.tsv", FORWARD_WIDE)
    lkp = _lookup_csv(tmp_path, lookup)
    outp = str(tmp_path / "out.tsv")
    assert main(["convert", inp, "--lookup", lkp, "--to", "AB", "-o", outp]) == 0
    lines = open(outp).read().splitlines()
    assert lines[1] == "snp1\tAA\tAB\tBB"


def test_cli_convert_unknown_snp_fails(tmp_path, lookup):
    inp = _write(tmp_path, "in.tsv", "SNP\ts1\nsnp_NOPE\tAA\n")
    lkp = _lookup_csv(tmp_path, lookup)
    outp = str(tmp_path / "out.tsv")
    assert main(["convert", inp, "--lookup", lkp, "--from", "AB",
                 "--to", "FORWARD", "-o", outp]) == 2


def test_console_script_entrypoint(tmp_path, lookup):
    # exercises the installed `soysnp-format` script via pip install -e
    inp = _write(tmp_path, "in.tsv", AB_WIDE)
    lkp = _lookup_csv(tmp_path, lookup)
    proc = subprocess.run(
        ["soysnp-format", "detect", inp, "--lookup", lkp],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Likely format: AB" in proc.stdout
