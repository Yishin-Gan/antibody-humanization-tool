"""
Flask app: Antibody Humanization Advisor.

Routes
------
GET  /                       Input form
POST /api/detect-vl-chain    Returns {chain_type, ok} for a pasted VL
POST /api/validate-germline  Returns {found, resolved, note}
POST /run                    Runs the pipeline (blocking, 30s-2min); redirects to /report/<id>
GET  /report/<id>            Interactive HTML report
GET  /api/report/<id>/data   JSON payload (used by tweak re-render path)
POST /api/report/<id>/rescore Re-score affected chain after residue edits
GET  /api/report/<id>/xlsx   Download structured workbook (Glossary + Sequences + Metrics)
GET  /api/report/<id>/fasta  Download all active mode sequences as FASTA
GET  /api/report/<id>/html   Download a self-contained copy of the report
"""

from __future__ import annotations

import json
import sys
from threading import Lock

sys.path.insert(0, "/workspace/antibody-humanization-tool")

from flask import Flask, render_template, request, redirect, url_for, jsonify, abort, Response

from web import runner, report_data, rescore, download


app = Flask(
    __name__,
    template_folder="templates",
    static_folder="static",
)
app.config["JSON_SORT_KEYS"] = False

# Single-user local app: stash run results in memory.
JOBS: dict = {}
JOBS_LOCK = Lock()


# ── input form ────────────────────────────────────────────────────────────────


@app.route("/")
def index():
    return render_template("input.html")


@app.post("/api/detect-vl-chain")
def api_detect_vl_chain():
    seq = request.json.get("seq", "") if request.is_json else request.form.get("seq", "")
    return jsonify(runner.detect_vl_chain_type(seq))


@app.post("/api/validate-germline")
def api_validate_germline():
    data = request.get_json(silent=True) or request.form
    name = (data.get("name") or "").strip()
    chain = (data.get("chain") or "H").strip().upper()
    if chain not in ("H", "K", "L"):
        return jsonify({"found": False, "resolved": None, "note": "invalid chain"})
    return jsonify(runner.validate_germline_name(name, chain))


# ── run pipeline ──────────────────────────────────────────────────────────────


@app.post("/run")
def run_pipeline():
    # Accept both classic form POST and JSON (from fetch-based submission).
    if request.is_json:
        body = request.get_json()
        f = body
        getlist = lambda k: body.get(k, []) if isinstance(body.get(k), list) else (
            [body[k]] if k in body and body[k] else [])
    else:
        f = request.form
        getlist = f.getlist
    modes = getlist("modes")
    if not modes:
        if request.is_json:
            return jsonify({"error": "Select at least one processing mode."}), 400
        return "<h2>Error</h2><p>Select at least one processing mode.</p>", 400

    mouse_vh = f.get("mouse_vh", "")
    mouse_vl = f.get("mouse_vl", "")
    vl_chain_type = f.get("vl_chain_type", "").strip()
    if vl_chain_type == "auto":
        vl_chain_type = ""

    kwargs = dict(
        mouse_vh=mouse_vh,
        mouse_vl=mouse_vl,
        modes=modes,
        vl_chain_type=vl_chain_type or None,
        preferred_germ_vh=f.get("preferred_germ_vh", "").strip() or None,
        preferred_germ_vl=f.get("preferred_germ_vl", "").strip() or None,
        lab_hu_vh=f.get("lab_hu_vh", ""),
        lab_hu_vl=f.get("lab_hu_vl", ""),
        lab_final_vh=f.get("lab_final_vh", ""),
        lab_final_vl=f.get("lab_final_vl", ""),
        lab_germ_vh=f.get("lab_germ_vh", "").strip() or None,
        lab_germ_vl=f.get("lab_germ_vl", "").strip() or None,
        run_structure=f.get("structure") == "on",
    )

    result = runner.run(**kwargs)
    payload = report_data.build(result)

    with JOBS_LOCK:
        JOBS[result.job_id] = {"result": result, "payload": payload}
    if request.is_json:
        return jsonify({"job_id": result.job_id, "report_url": url_for("show_report", job_id=result.job_id)})
    return redirect(url_for("show_report", job_id=result.job_id))


# ── report ────────────────────────────────────────────────────────────────────


def _get_job(job_id):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    if not j:
        abort(404, f"job {job_id} not found")
    return j


@app.get("/report/<job_id>")
def show_report(job_id):
    j = _get_job(job_id)
    payload_json = json.dumps(j["payload"], default=str)
    return render_template("report.html", payload_json=payload_json, job_id=job_id)


@app.get("/api/report/<job_id>/data")
def api_report_data(job_id):
    j = _get_job(job_id)
    return jsonify(j["payload"])


@app.post("/api/report/<job_id>/rescore")
def api_report_rescore(job_id):
    """Re-score one mode after residue-level edits.

    Body:
      {
        "mode":    "pipeline" | "preferred" | "sapiens",
        "edits":   {"VH": {"<linear_1based>": "<aa>", ...},
                    "VL": {"<linear_1based>": "<aa>", ...}},
        "structure": true|false   # optional, defaults to whether structure
                                  # was on in the original run
      }

    Returns:
      { "VH": {<rescored chain metrics>} | null,
        "VL": {<rescored chain metrics>} | null,
        "structure": {<conf_* fields>} | null }
    """
    j = _get_job(job_id)
    body = request.get_json(force=True)
    mode = body.get("mode")
    edits_all = body.get("edits") or {}
    payload = j["payload"]
    mode_data = payload["modes"].get(mode)
    if not mode_data:
        return jsonify({"error": f"mode {mode} not active"}), 400
    vl_chain_type = payload["chain_type"] or "K"

    run_structure = body.get("structure")
    if run_structure is None:
        run_structure = bool(payload.get("structure"))

    out = {"VH": None, "VL": None, "structure": None}
    edited_seqs = {}  # chain → final sequence (edited or original)

    for chain, ctype, seq_key, graft_key, germ_key, mouse_seq in [
        ("VH", "H",            "vh", "grafted_vh", "vh_germline", payload["mouse_vh"]),
        ("VL", vl_chain_type,  "vl", "grafted_vl", "vl_germline", payload["mouse_vl"]),
    ]:
        original_seq = mode_data[seq_key] or ""
        chain_edits = (edits_all.get(chain) or {})
        edited = rescore.apply_edits(original_seq, chain_edits)
        edited_seqs[chain] = edited
        if chain_edits:
            out[chain] = rescore.rescore_chain(
                edited_seq=edited,
                mouse_seq=mouse_seq,
                chain_type=ctype,
                germline_name=mode_data.get(germ_key) or "",
                grafted_seq=mode_data[graft_key] or "",
            )
            out[chain]["chain"] = chain

    # Structure prediction runs on the PAIRED Fv, so it sees whichever chain
    # was edited plus the current state of the other chain (which may also
    # have edits applied if both chains were edited in this call).
    if run_structure and edited_seqs["VH"] and edited_seqs["VL"]:
        out["structure"] = rescore.rescore_structure(
            edited_seqs["VH"], edited_seqs["VL"], clone_id=f"tweak_{job_id}")

    out["mode"] = mode
    out["structure_requested"] = bool(run_structure)
    return jsonify(out)


# ── downloads ─────────────────────────────────────────────────────────────────


@app.get("/api/report/<job_id>/xlsx")
def api_report_xlsx(job_id):
    j = _get_job(job_id)
    blob = download.build_xlsx(j["result"])
    return Response(
        blob,
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": f"attachment; filename=humanization_{job_id}.xlsx"})


@app.post("/api/report/<job_id>/tweak/fasta")
def api_tweak_fasta(job_id):
    """Apply tweak edits to one mode's sequences and return as FASTA."""
    j = _get_job(job_id)
    body = request.get_json(force=True)
    mode = body.get("mode")
    edits_all = body.get("edits") or {}
    payload = j["payload"]
    mode_data = payload["modes"].get(mode)
    if not mode_data:
        return jsonify({"error": f"mode {mode} not active"}), 400

    vh_edits = edits_all.get("VH") or {}
    vl_edits = edits_all.get("VL") or {}
    vh_seq = rescore.apply_edits(mode_data["vh"] or "", vh_edits)
    vl_seq = rescore.apply_edits(mode_data["vl"] or "", vl_edits)

    lines = []
    if vh_seq:
        lines += [f">{mode}_VH_{job_id}_tweaked", vh_seq]
    if vl_seq:
        lines += [f">{mode}_VL_{job_id}_tweaked", vl_seq]
    return Response("\n".join(lines) + "\n", mimetype="text/plain",
                    headers={"Content-Disposition": f"attachment; filename=tweaked_{mode}_{job_id}.fasta"})


@app.post("/api/report/<job_id>/tweak/xlsx")
def api_tweak_xlsx(job_id):
    """Apply edits + re-score the tweaked sequences and return an XLSX bundle.

    Structure prediction runs if it was enabled on the original submission.
    """
    j = _get_job(job_id)
    body = request.get_json(force=True)
    mode = body.get("mode")
    edits_all = body.get("edits") or {}
    payload = j["payload"]
    mode_data = payload["modes"].get(mode)
    if not mode_data:
        return jsonify({"error": f"mode {mode} not active"}), 400

    blob = download.build_tweak_xlsx(
        result=j["result"], payload=payload, mode=mode,
        vh_edits=edits_all.get("VH") or {},
        vl_edits=edits_all.get("VL") or {},
        run_structure=bool(payload.get("structure")),
    )
    return Response(
        blob,
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": f"attachment; filename=tweaked_{mode}_{job_id}.xlsx"})


@app.get("/api/report/<job_id>/fasta")
def api_report_fasta(job_id):
    j = _get_job(job_id)
    payload = j["payload"]
    lines = []
    for mode, d in payload["modes"].items():
        for chain in ("vh", "vl"):
            seq = d.get(chain) or ""
            if seq:
                lines.append(f">{mode}_{chain.upper()}_{payload['job_id']}")
                lines.append(seq)
    mouse_vh = payload.get("mouse_vh", "")
    mouse_vl = payload.get("mouse_vl", "")
    if mouse_vh:
        lines += [f">mouse_VH_{payload['job_id']}", mouse_vh]
    if mouse_vl:
        lines += [f">mouse_VL_{payload['job_id']}", mouse_vl]
    return Response("\n".join(lines) + "\n", mimetype="text/plain",
                    headers={"Content-Disposition": f"attachment; filename=sequences_{job_id}.fasta"})


@app.get("/api/report/<job_id>/html")
def api_report_html(job_id):
    j = _get_job(job_id)
    payload_json = json.dumps(j["payload"], default=str)
    rendered = render_template("report.html", payload_json=payload_json, job_id=job_id, standalone=True)
    return Response(rendered, mimetype="text/html",
                    headers={"Content-Disposition": f"attachment; filename=report_{job_id}.html"})


# ── dev launcher ──────────────────────────────────────────────────────────────


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False, threaded=False)
