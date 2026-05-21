// Antibody Humanization Advisor — interactive report

(() => {
const D = window.PAYLOAD;          // canonical payload from server
const JOB = window.JOB_ID;

const VERN = { VH: new Set(D.vernier_vh), VL: new Set(D.vernier_vl) };

const MODE_COLOR = {
  pipeline:  "var(--pipeline)",
  preferred: "var(--preferred)",
  sapiens:   "var(--sapiens)",
  lab:       "var(--lab)",
  mouse:     "var(--mouse)",
};
const MODE_LABEL_FALLBACK = {
  pipeline:  "Pipeline",
  preferred: "Preferred",
  sapiens:   "Sapiens",
  lab:       "Lab ref",
  mouse:     "Mouse",
};

// active mode → spec-aligned key list (sapiens has no germline identity etc.)
const ACTIVE = D.active_modes || [];

// Per-mode local edit state: { mode: { VH: {linearIdx1: aa}, VL: {...} } }
// linear index = 1-based position in the displayed sequence (not IMGT)
const EDITS = {};
ACTIVE.forEach(m => { EDITS[m] = { VH: {}, VL: {} }; });

// Per-mode local re-scored cache: { mode: { VH: data, VL: data } }
const RESCORED = {};

// =====================================================================
// helpers
// =====================================================================
function el(tag, attrs, ...children) {
  const e = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v;
    else if (k === "style") e.style.cssText = v;
    else if (k === "html") e.innerHTML = v;
    else if (k.startsWith("on") && typeof v === "function") e.addEventListener(k.slice(2), v);
    else if (k === "dataset") for (const [dk, dv] of Object.entries(v)) e.dataset[dk] = dv;
    else if (v != null && v !== false) e.setAttribute(k, v);
  }
  for (const c of children.flat()) {
    if (c == null || c === false) continue;
    if (typeof c === "string" || typeof c === "number") e.appendChild(document.createTextNode(c));
    else e.appendChild(c);
  }
  return e;
}
function pct(v, d=1) { return v == null ? "—" : (100 * v).toFixed(d) + "%"; }
function num(v, d=2) { return v == null ? "—" : (+v).toFixed(d); }
function intval(v) { return (v == null ? "—" : v); }
function label(mode) { return (D.modes[mode] && D.modes[mode].label) || MODE_LABEL_FALLBACK[mode]; }
function sortedActive() {
  const ord = ["pipeline", "preferred", "sapiens", "lab"];
  return ord.filter(m => ACTIVE.includes(m));
}

// Non-tweak tabs always render ORIGINAL data. The Tweak tab is the only
// place where re-scored values are displayed, so the user can flip back and
// forth (Summary/Sequence View/etc. = original; Tweak = current tweak).
function scoresFor(mode)  { return D.modes[mode].scores; }
function seqFor(mode, chain) {
  const key = chain === "VH" ? "vh" : "vl";
  return D.modes[mode][key] || "";
}
function regionsFor(mode, chain) {
  return chain === "VH" ? D.modes[mode].vh_regions : D.modes[mode].vl_regions;
}
function ppMapFor(mode, chain) {
  return chain === "VH" ? D.modes[mode].vh_pp : D.modes[mode].vl_pp;
}
function imgtMapFor(mode, chain) {
  if (mode === "mouse") return chain === "VH" ? D.mouse.vh_imgt : D.mouse.vl_imgt;
  return chain === "VH" ? D.modes[mode].vh_imgt : D.modes[mode].vl_imgt;
}
function mutableFor(mode, chain) {
  return chain === "VH" ? D.modes[mode].mutable_vh : D.modes[mode].mutable_vl;
}

// Tweak-side accessors — only renderTweak() uses these.
function tweakScoresFor(mode) {
  return (RESCORED[mode] && RESCORED[mode]._all) || null;
}
function tweakSeqFor(mode, chain) {
  if (RESCORED[mode] && RESCORED[mode][chain]) return RESCORED[mode][chain].seq;
  // No server re-score yet for this chain — apply edits to the original string
  const orig = seqFor(mode, chain);
  const ed = (EDITS[mode] || {})[chain] || {};
  if (!Object.keys(ed).length) return orig;
  const chars = orig.split("");
  for (const [k, v] of Object.entries(ed)) {
    const i = +k - 1;
    if (i >= 0 && i < chars.length) chars[i] = v;
  }
  return chars.join("");
}
function germDbFor(mode, chain) {
  return chain === "VH" ? D.modes[mode].vh_germ_db : D.modes[mode].vl_germ_db;
}
function germSeqFor(mode, chain) {
  return chain === "VH" ? D.modes[mode].germ_seq_vh : D.modes[mode].germ_seq_vl;
}

// =====================================================================
// sidebar legend
// =====================================================================
function renderLegend() {
  const legend = document.getElementById("legend");
  legend.innerHTML = "";
  const items = [
    ["mouse",     "Mouse (input)"],
    ["pipeline",  "Pipeline"],
    ["preferred", "Preferred"],
    ["sapiens",   "Sapiens"],
    ["lab",       "Lab ref"],
  ];
  for (const [mode, lbl] of items) {
    if (mode !== "mouse" && !ACTIVE.includes(mode)) continue;
    legend.appendChild(el("div", {class: "row"},
      el("span", {class: "dot", style: `background:${MODE_COLOR[mode]}`}),
      el("span", null, lbl)
    ));
  }
}

// =====================================================================
// Summary tab
// =====================================================================
function renderSummary() {
  const root = document.getElementById("tab-summary");
  root.innerHTML = "";
  if (D.error) {
    root.appendChild(el("div", {class: "section-card"},
      el("h2", null, "Pipeline error"),
      el("pre", {class: "error-banner"}, D.error)
    ));
    return;
  }

  // ── Header ───
  const modesLine = sortedActive().map(m => {
    if (m === "preferred") return `Preferred (${D.preferred_germ_vh || "?"} / ${D.preferred_germ_vl || "?"})`;
    if (m === "lab") return `Lab ref (${D.lab_germ_vh || "?"} / ${D.lab_germ_vl || "?"})`;
    if (m === "pipeline") return `Pipeline (${D.pipeline_germ_vh || "?"} / ${D.pipeline_germ_vl || "?"})`;
    return label(m);
  }).join(" · ");

  const summaryCard = el("div", {class: "section-card"},
    el("h2", null, "Summary"),
    el("p", {class: "muted"},
      `Mouse input: ${D.mouse_vh.length} aa VH · ${D.mouse_vl.length} aa VL (${
        D.chain_type === "K" ? "Kappa" : (D.chain_type === "L" ? "Lambda" : "?")
      })`),
    el("p", {class: "muted"}, "Modes run: " + modesLine),
  );
  root.appendChild(summaryCard);

  // ── Humanness bars ───
  const humCard = el("div", {class: "section-card"});
  humCard.appendChild(el("h2", null,
    el("span", {dataTip: "Fraction of 9-mer windows in framework regions that appear in ≥10% of human OAS antibodies. Framework only — CDRs excluded. Primary humanness benchmark."}, "Humanness — OASis FR")));
  humCard.appendChild(el("p", {class: "smol"}, "Framework-only OASis 9-mer identity (primary metric)"));
  const grid = el("div", {class: "oasis-bars-grid"});
  for (const chain of ["VH", "VL"]) {
    const block = el("div", {class: "bars-block"}, el("h3", null, chain));
    const all = sortedActive().concat(["mouse"]);
    for (const mode of all) {
      let val;
      let lbl;
      if (mode === "mouse") {
        val = chain === "VH" ? D.mouse.scores.oa_fr_vh : D.mouse.scores.oa_fr_vl;
        lbl = "Mouse";
      } else {
        const sc = scoresFor(mode);
        val = chain === "VH" ? sc.oa_fr_vh : sc.oa_fr_vl;
        lbl = label(mode);
      }
      const pctv = val != null ? (val * 100) : 0;
      block.appendChild(el("div", {class: "bar-row"},
        el("div", {class: "label"}, lbl),
        el("div", {class: `bar ${mode}`}, el("span", {style: `width:${Math.max(0,Math.min(100,pctv))}%`})),
        el("div", {class: "val"}, pct(val))
      ));
    }
    grid.appendChild(block);
  }
  humCard.appendChild(grid);
  root.appendChild(humCard);

  // ── Key metrics table ───
  const km = el("div", {class: "section-card"});
  km.appendChild(el("h2", null, "Key metrics"));
  const tbl = el("table", {class: "metrics"});
  // header
  const thead = el("thead");
  const hr = el("tr", null, el("th", null, "Metric"));
  for (const m of sortedActive()) hr.appendChild(el("th", {class: "num " + "swatch-" + m}, label(m)));
  thead.appendChild(hr);
  tbl.appendChild(thead);
  const tbody = el("tbody");

  function row(title, getter, fmt=pct, tipText=null) {
    const tr = el("tr");
    const labelCell = el("td", null, tipText ? el("span", {dataTip: tipText}, title) : title);
    tr.appendChild(labelCell);
    for (const m of sortedActive()) {
      const v = getter(m);
      tr.appendChild(el("td", {class: "num"}, v == null ? "n/a" : fmt(v)));
    }
    return tr;
  }
  function groupRow(text) {
    return el("tr", {class: "group-row"}, el("td", {colspan: 1 + sortedActive().length}, "─── " + text + " ───"));
  }

  tbody.appendChild(row("OASis FR VH", m => scoresFor(m).oa_fr_vh, pct,
    "OASis FR identity — fraction of FR 9-mers in ≥10% of human repertoire."));
  tbody.appendChild(row("OASis FR VL", m => scoresFor(m).oa_fr_vl, pct));
  tbody.appendChild(row("Germline FR identity VH", m => m === "sapiens" ? null : scoresFor(m).g_vh, pct,
    "Fraction of FR IMGT positions matching the canonical germline sequence. Higher = closer to human scaffold."));
  tbody.appendChild(row("Germline FR identity VL", m => m === "sapiens" ? null : scoresFor(m).g_vl, pct));

  tbody.appendChild(groupRow("VH Vernier"));
  tbody.appendChild(row("Vernier mutable VH",     m => m === "sapiens" ? null : scoresFor(m).vh_vern_mut,  intval,
    "Vernier FR positions where mouse ≠ germline — each is a mandatory decision point."));
  tbody.appendChild(row("Vernier back-mutated VH", m => m === "sapiens" ? null : scoresFor(m).vh_vern_back, intval));
  tbody.appendChild(row("Vernier humanized VH",    m => m === "sapiens" ? null : scoresFor(m).vh_vern_hum,  intval));
  tbody.appendChild(row("Vernier other VH",        m => m === "sapiens" ? null : scoresFor(m).vh_vern_other, intval,
    "Sapiens introduced a third residue (neither mouse nor germline) — unvalidated. High risk."));

  tbody.appendChild(groupRow("VL Vernier"));
  tbody.appendChild(row("Vernier mutable VL",      m => m === "sapiens" ? null : scoresFor(m).vl_vern_mut,  intval));
  tbody.appendChild(row("Vernier back-mutated VL", m => m === "sapiens" ? null : scoresFor(m).vl_vern_back, intval));
  tbody.appendChild(row("Vernier humanized VL",    m => m === "sapiens" ? null : scoresFor(m).vl_vern_hum,  intval));
  tbody.appendChild(row("Vernier other VL",        m => m === "sapiens" ? null : scoresFor(m).vl_vern_other, intval));

  tbody.appendChild(groupRow("All FR positions"));
  tbody.appendChild(row("Total mutable VH+VL",
    m => {
      const s = scoresFor(m);
      if (m === "sapiens") return null;
      return ((s.vh_fr_mut||0) + (s.vl_fr_mut||0));
    }, intval));
  tbody.appendChild(row("Total back-mutated VH+VL",
    m => {
      const s = scoresFor(m);
      if (m === "sapiens") return null;
      return ((s.vh_fr_back||0) + (s.vl_fr_back||0));
    }, intval));
  tbody.appendChild(row("Total humanized VH+VL",
    m => {
      const s = scoresFor(m);
      if (m === "sapiens") return null;
      return ((s.vh_fr_hum||0) + (s.vl_fr_hum||0));
    }, intval));

  tbody.appendChild(groupRow("Sapiens drift"));
  tbody.appendChild(row("Drift VH",
    m => {
      const s = scoresFor(m);
      if (m === "sapiens" || m === "lab") return null;
      if (!s.drift_vh_post) return "—";
      return s.drift_vh_flag
        ? `DRIFTED → ${s.drift_vh_post}`
        : `None (${s.drift_vh_post})`;
    }, v => v));
  tbody.appendChild(row("Drift VL",
    m => {
      const s = scoresFor(m);
      if (m === "sapiens" || m === "lab") return null;
      if (!s.drift_vl_post) return "—";
      return s.drift_vl_flag
        ? `DRIFTED → ${s.drift_vl_post}`
        : `None (${s.drift_vl_post})`;
    }, v => v));

  tbody.appendChild(groupRow("Solubility"));
  tbody.appendChild(row("CamSol FR VH", m => scoresFor(m).cs_fr_vh, v => num(v, 2),
    "Mean intrinsic CamSol solubility over FR residues. Higher = more soluble. Below 0 = aggregation-prone."));
  tbody.appendChild(row("CamSol FR VL", m => scoresFor(m).cs_fr_vl, v => num(v, 2)));

  tbl.appendChild(tbody);
  km.appendChild(tbl);
  root.appendChild(km);
}

// =====================================================================
// Sequence View tab
// =====================================================================
function buildRuler(length, imgtMap) {
  // One inline cell per residue, same 13px width as .aa. Cells are labelled
  // by their IMGT position so the Sequence View ruler agrees with the
  // Germline Identity alignment, Mutable Positions, Vernier set, etc. —
  // all of which speak IMGT. Tick every 10th IMGT position, plus the very
  // first cell.
  const row = el("div", {class: "seq-ruler-row"});
  for (let i = 0; i < length; i++) {
    const pos = imgtMap && imgtMap[i] != null ? imgtMap[i] : (i + 1);
    const isTick = (pos % 10 === 0) || i === 0;
    const cls = "ruler-cell" + (isTick ? " tick" : "");
    row.appendChild(el("span", {class: cls}, isTick ? String(pos) : ""));
  }
  return row;
}

function renderSequenceView() {
  const root = document.getElementById("tab-sequence");
  root.innerHTML = "";
  if (D.error) {
    root.appendChild(el("div", {class: "section-card error-banner"}, D.error));
    return;
  }

  for (const chain of ["VH", "VL"]) {
    const card = el("div", {class: "section-card"});
    card.appendChild(el("h2", null, chain));

    const block = el("div", {class: "seq-block"});

    // Pick the longest active sequence as the visual reference; its IMGT
    // mapping drives the ruler so labels match what the rest of the report
    // (alignment, mutable positions) speaks.
    let refSeq = "";
    let refRegions = null;
    let refImgt = null;
    for (const m of sortedActive()) {
      const s = seqFor(m, chain);
      if (s && s.length > refSeq.length) {
        refSeq = s;
        refRegions = regionsFor(m, chain);
        refImgt = imgtMapFor(m, chain);
      }
    }
    if (!refSeq) {
      refSeq = chain === "VH" ? D.mouse_vh : D.mouse_vl;
      refImgt = imgtMapFor("mouse", chain);
    }

    const rows = el("div", {class: "seq-rows"});
    // Ruler row first — same grid so columns align under the residue rows.
    rows.appendChild(el("div", {class: "row-label"}, ""));
    rows.appendChild(buildRuler(refSeq.length, refImgt));
    // Each active mode + mouse
    const seqList = sortedActive().concat(["mouse"]);
    // Reference for diff = pipeline if available, else first active
    const referenceForDiff = {};
    for (const m of sortedActive()) {
      const gdb = germDbFor(m, chain);
      referenceForDiff[m] = gdb;  // diff vs germline DB residue set for that mode
    }
    const pipelineGerm = (germDbFor("pipeline", chain) || germDbFor(sortedActive()[0], chain) || {});

    for (const m of seqList) {
      let seq = "";
      let regions = null;
      if (m === "mouse") {
        seq = chain === "VH" ? D.mouse_vh : D.mouse_vl;
        regions = chain === "VH" ? D.mouse.vh_regions : D.mouse.vl_regions;
      } else {
        seq = seqFor(m, chain);
        regions = regionsFor(m, chain);
      }
      if (!seq) continue;

      // Build per-position class info
      // - cdr: from regions
      // - vernier: from VERN[chain] + IMGT position (from server-emitted vh_imgt/vl_imgt)
      // - diff: residue differs from corresponding germline (for that mode)
      // - edited: in EDITS[m][chain]
      const inCDR = new Array(seq.length).fill(false);
      if (regions) for (const r of regions) if (r.t !== "FR") for (let i = r.start; i < r.end; i++) inCDR[i] = true;

      // linear → IMGT comes from the server (ANARCI-numbered)
      const posList = imgtMapFor(m, chain) || [];

      const rowDiv = el("div", {class: "seq-row"});
      rowDiv.dataset.mode = m;
      rowDiv.dataset.chain = chain;
      for (let i = 0; i < seq.length; i++) {
        const aa = seq[i];
        const imgtPos = posList[i] != null ? posList[i] : null;
        const isV = imgtPos != null && VERN[chain].has(imgtPos);
        const ed = EDITS[m] && EDITS[m][chain] && EDITS[m][chain][i + 1];
        const cls = ["aa"];
        if (inCDR[i]) cls.push("cdr");
        if (isV) cls.push("vernier");
        if (ed) cls.push("edited");
        if (m === "mouse") {
          if (imgtPos != null) {
            const ref = pipelineGerm[String(imgtPos)];
            if (ref && ref.aa && ref.aa !== aa) cls.push("diff", "mouse");
          }
        } else {
          const gdb = referenceForDiff[m] || {};
          if (imgtPos != null) {
            const ref = gdb[String(imgtPos)];
            if (ref && ref.aa && ref.aa !== aa) cls.push("diff", m);
          }
        }
        const span = el("span", {class: cls.join(" "), dataset: {imgt: imgtPos != null ? imgtPos : "", lin: i + 1}}, aa);
        if (window._TWEAK_ON && (m === "pipeline" || m === "preferred")) {
          span.classList.add("tweak-clickable");
          span.addEventListener("click", (ev) => openEditPopup(ev, m, chain, i + 1, imgtPos, aa));
        }
        rowDiv.appendChild(span);
      }
      rows.appendChild(el("div", {class: "row-label swatch-" + m}, label(m)));
      rows.appendChild(rowDiv);
    }
    block.appendChild(rows);
    card.appendChild(block);
    root.appendChild(card);
  }

  // Vernier marker tooltip
  document.querySelectorAll("#tab-sequence .aa.vernier").forEach(s => {
    if (!s.title) s.title = "Vernier zone position — structural support for CDR loops.";
  });
}

// =====================================================================
// Germline Identity tab
// =====================================================================
function renderGermline() {
  const root = document.getElementById("tab-germline");
  root.innerHTML = "";
  if (D.error) {
    root.appendChild(el("div", {class: "section-card error-banner"}, D.error));
    return;
  }

  // Identity summary table — one row per mode, VH and VL side-by-side
  const card = el("div", {class: "section-card"});
  card.appendChild(el("h2", null, "Germline identity"));
  const tbl = el("table", {class: "metrics fit"});
  tbl.appendChild(el("thead", null,
    el("tr", null,
      el("th", null, "Mode"),
      el("th", null, "VH germline"),
      el("th", {class: "num"}, "VH FR id"),
      el("th", null, "VL germline"),
      el("th", {class: "num"}, "VL FR id"),
    )));
  const tbody = el("tbody");
  for (const m of sortedActive()) {
    const d = D.modes[m];
    const sc = scoresFor(m);
    const ghv = d.vh_germline;
    const glv = d.vl_germline;
    const ghOk = ghv && (Object.keys(d.germ_seq_vh || {}).length > 0);
    const glOk = glv && (Object.keys(d.germ_seq_vl || {}).length > 0);
    tbody.appendChild(el("tr", null,
      el("td", {class: "swatch-" + m}, label(m)),
      el("td", null, ghv || el("span", {class: "smol"}, "—")),
      el("td", {class: "num"}, sc.g_vh != null ? pct(sc.g_vh) : (ghOk ? "—" : el("span", {class: "smol"}, "n/a"))),
      el("td", null, glv || el("span", {class: "smol"}, "—")),
      el("td", {class: "num"}, sc.g_vl != null ? pct(sc.g_vl) : (glOk ? "—" : el("span", {class: "smol"}, "n/a"))),
    ));
  }
  tbl.appendChild(tbody);
  card.appendChild(tbl);

  // Post-Sapiens drift summary
  for (const m of sortedActive()) {
    if (m === "sapiens" || m === "lab") continue;
    const sc = scoresFor(m);
    if (sc.drift_vh_post || sc.drift_vl_post) {
      const driftCard = el("div", {style: "margin-top:0.6em"});
      driftCard.appendChild(el("h3", null, `${label(m)} — post-Sapiens drift`));
      if (sc.drift_vh_post) {
        const cls = sc.drift_vh_flag ? "drift-line warn" : "drift-line";
        const txt = sc.drift_vh_flag
          ? `Post-Sapiens VH → ${sc.drift_vh_post}  Δ=${sc.drift_vh_delta != null ? (sc.drift_vh_delta >= 0 ? '+' : '') + sc.drift_vh_delta.toFixed(3) : "?"}  ⚠ DRIFTED`
          : `Post-Sapiens VH → ${sc.drift_vh_post}  Δ=${sc.drift_vh_delta != null ? (sc.drift_vh_delta >= 0 ? '+' : '') + sc.drift_vh_delta.toFixed(3) : "?"}  (no drift)`;
        driftCard.appendChild(el("div", {class: cls}, txt));
      }
      if (sc.drift_vl_post) {
        const cls = sc.drift_vl_flag ? "drift-line warn" : "drift-line";
        const txt = sc.drift_vl_flag
          ? `Post-Sapiens VL → ${sc.drift_vl_post}  Δ=${sc.drift_vl_delta != null ? (sc.drift_vl_delta >= 0 ? '+' : '') + sc.drift_vl_delta.toFixed(3) : "?"}  ⚠ DRIFTED`
          : `Post-Sapiens VL → ${sc.drift_vl_post}  Δ=${sc.drift_vl_delta != null ? (sc.drift_vl_delta >= 0 ? '+' : '') + sc.drift_vl_delta.toFixed(3) : "?"}  (no drift)`;
        driftCard.appendChild(el("div", {class: cls}, txt));
      }
      card.appendChild(driftCard);
    }
  }

  root.appendChild(card);

  // ── Per-mode FR alignment panels ───
  const alignCard = el("div", {class: "section-card"});
  alignCard.appendChild(el("h2", null, "Framework alignment"));
  const tabs = el("div", {class: "align-tabs"});
  const panes = el("div");
  let firstPaneId = null;
  for (const m of sortedActive()) {
    const paneId = "align-" + m;
    if (firstPaneId == null) firstPaneId = paneId;
    tabs.appendChild(el("button", {
      class: "align-tab-btn" + (m === sortedActive()[0] ? " active" : ""),
      dataset: {pane: paneId},
      onclick: (ev) => {
        document.querySelectorAll(".align-tab-btn").forEach(b => b.classList.remove("active"));
        ev.target.classList.add("active");
        document.querySelectorAll("#tab-germline .align-pane").forEach(p => p.hidden = true);
        document.getElementById(paneId).hidden = false;
      },
    }, label(m)));
    const pane = el("div", {class: "align-pane", id: paneId, hidden: m !== sortedActive()[0]});
    for (const chain of ["VH", "VL"]) {
      pane.appendChild(buildAlignmentBlock(m, chain));
    }
    panes.appendChild(pane);
  }
  alignCard.appendChild(tabs);
  alignCard.appendChild(panes);
  alignCard.appendChild(el("p", {class: "smol"},
    "FR identity is computed against the canonical germline allele in the database. The alignment shows the scaffold sequence used during grafting, which may use a specific allele differing slightly from the canonical entry. Minor discrepancies between the identity percentage and the visible diff count are expected and reflect allele-level variation."));
  root.appendChild(alignCard);
}

function buildAlignmentBlock(mode, chain) {
  const d = D.modes[mode];
  const sc = scoresFor(mode);
  const block = el("div", {class: "align-block"});
  const germName = chain === "VH" ? d.vh_germline : d.vl_germline;
  if (!germName) {
    block.appendChild(el("p", {class: "smol"}, `${chain}: no germline reference for this mode.`));
    return block;
  }
  const identity = chain === "VH" ? sc.g_vh : sc.g_vl;
  block.appendChild(el("p", null,
    el("span", {class: "swatch-" + mode}, label(mode)), `  vs  `,
    el("span", {class: "swatch-germline"}, germName),
    `  ·  FR identity: `,
    el("strong", null, identity != null ? pct(identity) : "—")
  ));

  const germRegions = germSeqFor(mode, chain);
  const modeSeq = chain === "VH" ? d.vh : d.vl;
  const posList = imgtMapFor(mode, chain) || [];

  // Group residues by IMGT position — NOT by find_regions's substring-based
  // linear boundaries. find_regions can mis-assign boundary residues (e.g.
  // an FR3-start residue that abnumber numbers as IMGT 66 gets stuck in the
  // CDR2 substring), which would drop the first residue of an FR row.
  // Walking by IMGT guarantees the alignment matches abnumber numbering.
  const FR_REGIONS = ["FR1", "FR2", "FR3", "FR4"];
  const buckets = {FR1:[], FR2:[], FR3:[], FR4:[]};
  for (let i = 0; i < modeSeq.length; i++) {
    const imgti = posList[i];
    if (imgti == null) continue;
    let key;
    if      (imgti <= 26)  key = "FR1";
    else if (imgti <= 38)  continue;       // CDR1
    else if (imgti <= 55)  key = "FR2";
    else if (imgti <= 65)  continue;       // CDR2
    else if (imgti <= 104) key = "FR3";
    else if (imgti <= 117) continue;       // CDR3
    else                   key = "FR4";
    const germ = germRegions[String(imgti)];
    buckets[key].push({lin: i, imgt: imgti, aa: modeSeq[i], germ: germ ? germ.aa : null});
  }

  // Friendly germline display name for the row label (e.g. "IGHV1-69*02")
  const germLabel = (chain === "VH" ? d.vh_germline : d.vl_germline) || "germline";

  for (const region of FR_REGIONS) {
    const items = buckets[region];
    if (!items.length) continue;
    const sub = el("div", {class: "align-region"});
    sub.appendChild(el("p", {class: "smol align-region-title"},
      `${region} (IMGT ${items[0].imgt}–${items[items.length-1].imgt}):`));

    // Ruler row: one cell per residue, tick every 10 IMGT positions
    const rulerRow = el("div", {class: "align-row align-ruler"},
      el("span", {class: "align-label"}, ""));
    const rulerCells = el("div", {class: "align-cells"});
    for (const it of items) {
      const isTick = (it.imgt % 10 === 0) || it === items[0];
      rulerCells.appendChild(el("span",
        {class: "align-cell ruler" + (isTick ? " tick" : "")},
        isTick ? String(it.imgt) : ""));
    }
    rulerRow.appendChild(rulerCells);
    sub.appendChild(rulerRow);

    // Mode row
    const modeRow = el("div", {class: "align-row"},
      el("span", {class: "align-label swatch-" + mode}, label(mode)));
    const modeCells = el("div", {class: "align-cells"});
    // Germline row — label uses the actual germline name (e.g. "IGHV1-69*02")
    const germRow = el("div", {class: "align-row"},
      el("span", {class: "align-label swatch-germline", title: germLabel}, germLabel));
    const germCells = el("div", {class: "align-cells"});
    // Diff row
    const diffRow = el("div", {class: "align-row diff-row"},
      el("span", {class: "align-label"}, "Diff"));
    const diffCells = el("div", {class: "align-cells"});

    const changes = [];
    for (const it of items) {
      const isV = VERN[chain].has(it.imgt);
      const modeAA = it.aa || "-";
      const germAA = it.germ != null ? it.germ : "-";
      const isDiff = it.germ != null && it.aa != null && it.aa !== it.germ;

      const modeCell = el("span", {class: "align-cell"}, modeAA);
      if (isDiff) modeCell.classList.add("diff", mode);
      if (isV) modeCell.classList.add("vernier");
      modeCells.appendChild(modeCell);

      // Germline cell: only colored if it differs from mode (per request).
      // Otherwise muted so the eye is drawn to the diffs.
      const germCell = el("span", {class: "align-cell"}, germAA);
      if (isDiff) germCell.classList.add("germ-diff");
      else        germCell.classList.add("germ-match");
      if (isV) germCell.classList.add("vernier");
      germCells.appendChild(germCell);

      diffCells.appendChild(el("span",
        {class: "align-cell diff-mark" + (isDiff ? " on" : "")},
        isDiff ? "^" : ""));

      if (isDiff) changes.push(`pos${it.imgt}(${germAA}→${modeAA})${isV?"(V)":""}`);
    }
    modeRow.appendChild(modeCells);
    germRow.appendChild(germCells);
    diffRow.appendChild(diffCells);
    sub.appendChild(modeRow);
    sub.appendChild(germRow);
    if (changes.length) sub.appendChild(diffRow);

    // Changes line with wider separation between tokens (spec request)
    const changesEl = el("div", {class: "align-changes"});
    if (changes.length === 0) {
      changesEl.appendChild(el("span", {class: "smol"}, "identical"));
    } else {
      changes.forEach((c, idx) => {
        if (idx > 0) changesEl.appendChild(el("span", {class: "change-sep"}, " · "));
        changesEl.appendChild(el("span", {class: "change-token"}, c));
      });
    }
    sub.appendChild(changesEl);
    block.appendChild(sub);
  }
  return block;
}

// =====================================================================
// Mutable Positions tab
// =====================================================================
// Build a residue cell with colour applied directly to the residue character,
// matching the existing report HTML's style: coloured letter + optional symbol
// suffix ('°' not-mutable, '!' Sapiens-changed-non-mutable, '*' Sapiens-kept-mouse).
const RES_SUFFIX = {
  not_mutable:                 "°",
  sapiens_changed_non_mutable: "!",
  sapiens_kept_mouse:          "*",
};
const RES_TIP = {
  back_mutated:                "back-mutated to mouse",
  humanized:                   "kept as germline",
  humanized_by_sapiens:        "humanised by Sapiens (no grafted baseline)",
  other_substitution:          "third residue — neither mouse nor germline",
  sapiens_changed_non_mutable: "Sapiens changed a non-mutable position",
  not_mutable:                 "mouse = germline — not mutable",
  sapiens_kept_mouse:          "Sapiens kept mouse residue",
};
function residueCell(aa, status) {
  const txt = (aa || "·") + (RES_SUFFIX[status] || "");
  return el("span",
    {class: "res res-" + (status || "empty"), title: RES_TIP[status] || ""},
    txt);
}

function legendSwatch(status, label, demo) {
  return el("span", {class: "legend-chip"},
    el("span", {class: "res res-" + status}, (demo || "A") + (RES_SUFFIX[status] || "")),
    el("span", null, label));
}

function renderMutable() {
  const root = document.getElementById("tab-mutable");
  root.innerHTML = "";
  if (D.error) {
    root.appendChild(el("div", {class: "section-card error-banner"}, D.error));
    return;
  }

  // Aggregate across modes — union of positions per chain
  // Build {chain: {pos: {pos, vernier, mouse, perMode: {mode: {aa, status}}}}}
  const agg = { VH: {}, VL: {} };
  for (const m of sortedActive()) {
    for (const chain of ["VH", "VL"]) {
      const rows = mutableFor(m, chain);
      for (const r of rows) {
        if (!agg[chain][r.pos]) {
          agg[chain][r.pos] = {
            pos: r.pos, vernier: r.vernier, mouse: r.mouse, perMode: {},
          };
        }
        agg[chain][r.pos].vernier = agg[chain][r.pos].vernier || r.vernier;
        agg[chain][r.pos].perMode[m] = {aa: r.query, status: r.status, grafted: r.grafted};
      }
    }
  }

  // Header card
  const head = el("div", {class: "section-card"});
  head.appendChild(el("h2", null, "Mutable positions"));
  head.appendChild(el("p", {class: "muted"},
    "All FR positions where mouse ≠ germline. Vernier-zone positions are higher structural risk."));

  // Totals row — fit-content so columns don't sprawl across full width
  const totalsTbl = el("table", {class: "metrics fit"});
  const totHead = el("tr", null, el("th", null, ""));
  for (const m of sortedActive()) totHead.appendChild(el("th", {class: "num swatch-" + m}, label(m)));
  totalsTbl.appendChild(el("thead", null, totHead));
  const totBody = el("tbody");
  function totRow(title, getter) {
    const tr = el("tr");
    tr.appendChild(el("td", null, title));
    for (const m of sortedActive()) {
      const v = getter(m);
      tr.appendChild(el("td", {class: "num"}, v == null ? "n/a" : v));
    }
    return tr;
  }
  totBody.appendChild(totRow("Total mutable VH+VL",
    m => {
      if (m === "sapiens") return null;
      const s = scoresFor(m);
      return (s.vh_fr_mut || 0) + (s.vl_fr_mut || 0);
    }));
  totBody.appendChild(totRow("Total back-mutated VH+VL",
    m => {
      if (m === "sapiens") return null;
      const s = scoresFor(m);
      return (s.vh_fr_back || 0) + (s.vl_fr_back || 0);
    }));
  totBody.appendChild(totRow("Total humanized VH+VL",
    m => {
      if (m === "sapiens") return null;
      const s = scoresFor(m);
      return (s.vh_fr_hum || 0) + (s.vl_fr_hum || 0);
    }));
  totalsTbl.appendChild(totBody);
  head.appendChild(totalsTbl);

  head.appendChild(el("div", {class: "mut-legend"},
    legendSwatch("back_mutated",               "back-mut"),
    legendSwatch("humanized",                  "humanized"),
    legendSwatch("other_substitution",         "third residue"),
    legendSwatch("sapiens_changed_non_mutable","Sapiens changed non-mutable"),
    legendSwatch("not_mutable",                "not mutable"),
    el("span", {class: "legend-chip"},
      el("span", {style: "display:inline-block;width:6px;height:14px;background:var(--vernier);border-radius:1px"}),
      el("span", null, "= Vernier")),
    el("span", {class: "legend-chip"},
      el("span", {style: "display:inline-block;width:18px;height:10px;background:rgba(245,165,35,0.25);border-radius:2px"}),
      el("span", null, "= lab/pipeline conflict")),
  ));
  root.appendChild(head);

  // 2×2 grid order per request:
  //   row1: VH Vernier  | VL Vernier
  //   row2: VH Non-Vern | VL Non-Vern
  const grid = el("div", {class: "mut-grid"});

  for (const kind of ["vern", "non"]) {
    for (const chain of ["VH", "VL"]) {
      const card = el("div", {class: "section-card"});
      card.appendChild(el("h2", null,
        `${chain} — ${kind === "vern" ? "Vernier zone" : "Non-Vernier FR"}`
      ));
      const positions = Object.values(agg[chain])
        .filter(r => (kind === "vern") === !!r.vernier)
        .sort((a, b) => a.pos - b.pos);

      if (positions.length === 0) {
        card.appendChild(el("div", {class: "mut-empty"}, "No positions"));
        grid.appendChild(card);
        continue;
      }

      const tbl = el("table", {class: "metrics"});
      const hr = el("tr", null,
        el("th", null, "Pos"),
        el("th", {class: "swatch-mouse"}, "Mouse"),
      );
      for (const m of sortedActive()) hr.appendChild(el("th", {class: "swatch-" + m}, label(m)));
      hr.appendChild(el("th", null, "Note"));
      tbl.appendChild(el("thead", null, hr));

      const body = el("tbody");
      for (const row of positions) {
        const tr = el("tr");
        if (row.vernier) tr.classList.add("vern");

        // Detect lab/pipeline conflict
        const pipeOrPref = row.perMode.pipeline || row.perMode.preferred;
        const labOne = row.perMode.lab;
        const labBack = labOne && labOne.status === "back_mutated";
        const pipeHum = pipeOrPref && pipeOrPref.status === "humanized";
        if (labBack && pipeHum) tr.classList.add("conflict");

        tr.appendChild(el("td", null, row.pos));
        tr.appendChild(el("td", null,
          el("span", {class: "res res-empty swatch-mouse"}, row.mouse)));
        for (const m of sortedActive()) {
          const cell = el("td");
          const v = row.perMode[m];
          if (v) {
            cell.appendChild(residueCell(v.aa, v.status));
          } else {
            cell.appendChild(el("span", {class: "smol"}, "—"));
          }
          tr.appendChild(cell);
        }
        // Notes
        const notes = [];
        if (labBack && row.perMode.pipeline && row.perMode.pipeline.status === "humanized")
          notes.push("Lab back-mutated; pipeline kept human");
        if (labBack && row.perMode.preferred && row.perMode.preferred.status === "humanized")
          notes.push("Lab back-mutated; preferred kept human");
        if (row.perMode.pipeline && row.perMode.preferred &&
            row.perMode.pipeline.aa !== row.perMode.preferred.aa)
          notes.push(`Germline choice affects this position: Pipe=${row.perMode.pipeline.aa} Pref=${row.perMode.preferred.aa}`);
        for (const m of sortedActive()) {
          const v = row.perMode[m];
          if (!v) continue;
          if (v.status === "sapiens_changed_non_mutable") notes.push(`${label(m)}: Sapiens changed non-mutable position`);
          if (v.status === "other_substitution") notes.push(`${label(m)}: third residue — neither mouse nor germline`);
        }
        tr.appendChild(el("td", {class: "smol"}, notes.length ? notes.join("; ") : "—"));
        body.appendChild(tr);
      }
      tbl.appendChild(body);
      card.appendChild(tbl);
      grid.appendChild(card);
    }
  }
  root.appendChild(grid);

  // Suggested back-mutations panel (lab vs pipeline conflicts)
  if (ACTIVE.includes("lab") && (ACTIVE.includes("pipeline") || ACTIVE.includes("preferred"))) {
    const referenceMode = ACTIVE.includes("pipeline") ? "pipeline" : "preferred";
    const conflicts = { VH: [], VL: [] };
    for (const chain of ["VH", "VL"]) {
      for (const row of Object.values(agg[chain])) {
        const lab = row.perMode.lab;
        const ref = row.perMode[referenceMode];
        if (lab && ref && lab.status === "back_mutated" && ref.status === "humanized") {
          conflicts[chain].push({pos: row.pos, mouse: row.mouse, ref_aa: ref.aa, lab_aa: lab.aa});
        }
      }
    }
    if (conflicts.VH.length || conflicts.VL.length) {
      const sg = el("div", {class: "suggest-box"});
      sg.appendChild(el("h3", null, "▼ Suggested back-mutations — where lab back-mutated but " + referenceMode + " kept human"));
      const lines = [];
      for (const chain of ["VH", "VL"]) {
        if (!conflicts[chain].length) continue;
        const items = conflicts[chain].map(c => `pos${c.pos} (${c.ref_aa}→${c.lab_aa})`).join("  ");
        sg.appendChild(el("p", null, el("strong", {class: "swatch-mouse"}, chain + ": "), items));
        const copyStr = conflicts[chain].map(c => `${c.pos}:${c.lab_aa}`).join(",");
        sg.appendChild(el("button", {class: "btn small", onclick: () => navigator.clipboard.writeText(copyStr)},
          `Copy ${chain} mutations`));
      }
      sg.appendChild(el("p", {class: "smol"}, "Applying these mutations to the " + referenceMode + " sequence would bring it closer to the lab's humanization approach."));
      root.appendChild(sg);
    }
  }
}

// =====================================================================
// Feature Metrics tab
// =====================================================================
function renderFeatures() {
  const root = document.getElementById("tab-features");
  root.innerHTML = "";
  if (D.error) {
    root.appendChild(el("div", {class: "section-card error-banner"}, D.error));
    return;
  }

  // OASis bars (4 columns × N modes)
  const card = el("div", {class: "section-card"});
  card.appendChild(el("h2", null, "OASis humanness — bars"));
  const grp = el("div", {class: "oasis-bars-grid"});
  for (const chain of ["VH", "VL"]) {
    const block = el("div", {class: "bars-block"});
    block.appendChild(el("h3", null, chain));
    function bg(label, getter) {
      const wrap = el("div", {style: "margin: 0.4em 0"});
      wrap.appendChild(el("div", {class: "smol"}, label));
      for (const m of sortedActive()) {
        const v = getter(m);
        const pctv = v != null ? (v * 100) : 0;
        wrap.appendChild(el("div", {class: "bar-row"},
          el("div", {class: "label"}, MODE_LABEL_FALLBACK[m] || m),
          el("div", {class: `bar ${m}`}, el("span", {style: `width:${Math.max(0,Math.min(100,pctv))}%`})),
          el("div", {class: "val"}, pct(v))
        ));
      }
      return wrap;
    }
    block.appendChild(bg("Overall",  m => chain === "VH" ? scoresFor(m).oa_vh   : scoresFor(m).oa_vl));
    block.appendChild(bg("FR only",  m => chain === "VH" ? scoresFor(m).oa_fr_vh : scoresFor(m).oa_fr_vl));
    block.appendChild(bg("CDR only", m => chain === "VH" ? scoresFor(m).oa_cdr_vh : scoresFor(m).oa_cdr_vl));
    block.appendChild(bg("Germline FR ID", m => m === "sapiens" ? null : (chain === "VH" ? scoresFor(m).g_vh : scoresFor(m).g_vl)));
    grp.appendChild(block);
  }
  card.appendChild(grp);
  root.appendChild(card);

  // Full metrics table
  const tcard = el("div", {class: "section-card"});
  tcard.appendChild(el("h2", null, "All metrics"));
  const tbl = el("table", {class: "metrics"});
  const head = el("tr", null, el("th", null, "Feature"));
  const modes = sortedActive();
  for (const m of modes) head.appendChild(el("th", {class: "num swatch-" + m}, label(m)));
  if (modes.length >= 2) head.appendChild(el("th", null, "Winner"));
  tbl.appendChild(el("thead", null, head));
  const body = el("tbody");

  function pushGroup(title) {
    body.appendChild(el("tr", {class: "group-row"},
      el("td", {colspan: head.children.length}, "─── " + title + " ───")));
  }
  function pushRow(rowTitle, getter, fmt=pct, direction="higher") {
    const tr = el("tr");
    tr.appendChild(el("td", null, rowTitle));
    const vals = modes.map(m => getter(m));
    for (const v of vals) tr.appendChild(el("td", {class: "num"}, v == null ? "n/a" : fmt(v)));
    if (modes.length >= 2) {
      let winnerIdx = -1;
      if (direction === "higher" || direction === "lower") {
        const cmp = direction === "higher" ? ((a, b) => a > b) : ((a, b) => a < b);
        let best = null;
        for (let i = 0; i < vals.length; i++) {
          if (vals[i] == null) continue;
          if (best == null || cmp(vals[i], best)) { best = vals[i]; winnerIdx = i; }
        }
        if (winnerIdx === -1) tr.appendChild(el("td", null, "—"));
        else tr.appendChild(el("td", {class: "swatch-" + modes[winnerIdx]}, "▲ " + label(modes[winnerIdx])));
      } else {
        tr.appendChild(el("td", null, "—"));
      }
    }
    return tr;
  }

  pushGroup("Humanness — OASis 9-mer prevalence (higher = more human)");
  body.appendChild(pushRow("OASis VH overall",  m => scoresFor(m).oa_vh,    pct));
  body.appendChild(pushRow("OASis FR VH",       m => scoresFor(m).oa_fr_vh, pct));
  body.appendChild(pushRow("OASis CDR VH",      m => scoresFor(m).oa_cdr_vh, pct));
  body.appendChild(pushRow("Germline FR id VH", m => m === "sapiens" ? null : scoresFor(m).g_vh, pct));
  body.appendChild(pushRow("OASis VL overall",  m => scoresFor(m).oa_vl,    pct));
  body.appendChild(pushRow("OASis FR VL",       m => scoresFor(m).oa_fr_vl, pct));
  body.appendChild(pushRow("OASis CDR VL",      m => scoresFor(m).oa_cdr_vl, pct));
  body.appendChild(pushRow("Germline FR id VL", m => m === "sapiens" ? null : scoresFor(m).g_vl, pct));

  const anyStructure = modes.some(m => scoresFor(m).cf_vh != null);
  if (anyStructure) {
    pushGroup("Structure — ABodyBuilder2 (higher = more confident)");
    body.appendChild(pushRow("Conf mean VH",  m => scoresFor(m).cf_vh,    v => num(v, 3)));
    body.appendChild(pushRow("Conf FR VH",    m => scoresFor(m).cf_fr_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR VH",   m => scoresFor(m).cf_cdr_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf mean VL",  m => scoresFor(m).cf_vl,    v => num(v, 3)));
    body.appendChild(pushRow("Conf FR VL",    m => scoresFor(m).cf_fr_vl, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR VL",   m => scoresFor(m).cf_cdr_vl, v => num(v, 3)));
    body.appendChild(pushRow("Conf min VH",   m => scoresFor(m).cf_min_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf min VL",   m => scoresFor(m).cf_min_vl, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR1 VH",  m => scoresFor(m).cf_cdr1_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR2 VH",  m => scoresFor(m).cf_cdr2_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR3 VH",  m => scoresFor(m).cf_cdr3_vh, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR1 VL",  m => scoresFor(m).cf_cdr1_vl, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR2 VL",  m => scoresFor(m).cf_cdr2_vl, v => num(v, 3)));
    body.appendChild(pushRow("Conf CDR3 VL",  m => scoresFor(m).cf_cdr3_vl, v => num(v, 3)));
  }

  pushGroup("Solubility — CamSol intrinsic (higher = more soluble)");
  body.appendChild(pushRow("CamSol VH",       m => scoresFor(m).cs_vh,     v => num(v, 3)));
  body.appendChild(pushRow("CamSol FR VH",    m => scoresFor(m).cs_fr_vh,  v => num(v, 3)));
  body.appendChild(pushRow("CamSol CDR VH",   m => scoresFor(m).cs_cdr_vh, v => num(v, 3)));
  body.appendChild(pushRow("CamSol VL",       m => scoresFor(m).cs_vl,     v => num(v, 3)));
  body.appendChild(pushRow("CamSol FR VL",    m => scoresFor(m).cs_fr_vl,  v => num(v, 3)));
  body.appendChild(pushRow("CamSol CDR VL",   m => scoresFor(m).cs_cdr_vl, v => num(v, 3)));
  body.appendChild(pushRow("Hotspots VH",     m => scoresFor(m).hs_vh,     intval, "lower"));
  body.appendChild(pushRow("Hotspots FR VH",  m => scoresFor(m).hs_fr_vh,  intval, "lower"));
  body.appendChild(pushRow("Hotspots CDR VH", m => scoresFor(m).hs_cdr_vh, intval, "neither"));
  body.appendChild(pushRow("Hotspots VL",     m => scoresFor(m).hs_vl,     intval, "lower"));
  body.appendChild(pushRow("Hotspots FR VL",  m => scoresFor(m).hs_fr_vl,  intval, "lower"));
  body.appendChild(pushRow("Hotspots CDR VL", m => scoresFor(m).hs_cdr_vl, intval, "neither"));

  pushGroup("Physicochemical");
  body.appendChild(pushRow("pI Fv",          m => scoresFor(m).pi, v => num(v, 2), "neither"));
  body.appendChild(pushRow("Net charge Fv",  m => scoresFor(m).ch, v => num(v, 2), "neither"));

  pushGroup("FR liabilities (lower = fewer risks)");
  body.appendChild(pushRow("Deamid FR VH",   m => scoresFor(m).lia.vh_df, intval, "lower"));
  body.appendChild(pushRow("Deamid FR VL",   m => scoresFor(m).lia.vl_df, intval, "lower"));
  body.appendChild(pushRow("Oxidation FR VH",m => scoresFor(m).lia.vh_of, intval, "lower"));
  body.appendChild(pushRow("Oxidation FR VL",m => scoresFor(m).lia.vl_of, intval, "lower"));
  body.appendChild(pushRow("Isomer FR VH",   m => scoresFor(m).lia.vh_if, intval, "lower"));
  body.appendChild(pushRow("Isomer FR VL",   m => scoresFor(m).lia.vl_if, intval, "lower"));

  pushGroup("CDR liabilities — intrinsic, identical across all pipeline sequences");
  body.appendChild(pushRow("Deamid CDR VH",   m => scoresFor(m).lia.vh_dc, intval, "neither"));
  body.appendChild(pushRow("Oxidation CDR VH",m => scoresFor(m).lia.vh_oc, intval, "neither"));
  body.appendChild(pushRow("Deamid CDR VL",   m => scoresFor(m).lia.vl_dc, intval, "neither"));
  body.appendChild(pushRow("Oxidation CDR VL",m => scoresFor(m).lia.vl_oc, intval, "neither"));

  tbl.appendChild(body);
  tcard.appendChild(tbl);
  root.appendChild(tcard);
}

// =====================================================================
// Tweak tab — current edits + side-by-side comparison vs original
// =====================================================================
function activeTweakMode() {
  // Pick the first editable mode that has edits or rescored results.
  const candidates = sortedActive().filter(m => m !== "sapiens" && m !== "lab");
  for (const m of candidates) {
    if (RESCORED[m] && RESCORED[m]._all) return m;
    if (EDITS[m] && (Object.keys(EDITS[m].VH).length || Object.keys(EDITS[m].VL).length)) return m;
  }
  return candidates[0] || null;
}

function fmtDelta(o, t, fmt) {
  if (o == null || t == null) return "—";
  const d = t - o;
  const sign = d > 0 ? "+" : "";
  if (Math.abs(d) < 1e-9) return el("span", {class: "delta-zero"}, "±0");
  const cls = d > 0 ? "delta-pos" : "delta-neg";
  const txt = sign + (fmt === pct ? (d*100).toFixed(2) + " pp" : (Math.abs(d) >= 1 ? d.toFixed(2) : d.toFixed(3)));
  return el("span", {class: cls}, txt);
}

function renderTweak() {
  const root = document.getElementById("tab-tweak");
  root.innerHTML = "";

  const mode = activeTweakMode();
  if (!mode) {
    root.appendChild(el("div", {class: "section-card"},
      el("p", {class: "muted"}, "No editable modes available.")));
    return;
  }

  const hasEdits = EDITS[mode] && (Object.keys(EDITS[mode].VH).length || Object.keys(EDITS[mode].VL).length);
  const tw = tweakScoresFor(mode);

  // ── Header card with mode picker + state + downloads ────────────────────
  const head = el("div", {class: "section-card"});
  head.appendChild(el("h2", null, "Tweak"));
  head.appendChild(el("p", {class: "muted"},
    "Edit residues in ", el("strong", null, "Sequence View"),
    " (Enable editing in the sidebar). Other tabs always show the original run; this tab shows the current tweak and Δ vs original."));

  const ctlRow = el("div", {style: "display:flex;align-items:center;gap:1em;flex-wrap:wrap;margin:.8em 0"});
  // Mode picker if more than one editable mode
  const editable = sortedActive().filter(m => m !== "sapiens" && m !== "lab");
  if (editable.length > 1) {
    const sel = el("select", {id: "tweak-mode-tab", style: "background:var(--field);color:var(--t1);border:1px solid var(--card-bd);border-radius:5px;padding:.3em .5em;font-size:.85em"});
    for (const m of editable) {
      const opt = el("option", {value: m}, label(m));
      if (m === mode) opt.selected = true;
      sel.appendChild(opt);
    }
    sel.addEventListener("change", () => {
      document.getElementById("tweak-mode").value = sel.value;
      renderTweak();
    });
    ctlRow.appendChild(el("span", {class: "smol"}, "Mode:"));
    ctlRow.appendChild(sel);
  }
  ctlRow.appendChild(el("span", {class: "swatch-" + mode, style: "font-weight:600"}, "● " + label(mode)));

  const editCount = (Object.keys((EDITS[mode]||{}).VH || {}).length) +
                    (Object.keys((EDITS[mode]||{}).VL || {}).length);
  ctlRow.appendChild(el("span", {class: "smol"},
    `${editCount} position${editCount === 1 ? "" : "s"} edited` +
    (tw ? " · re-scored" : (editCount ? " · not yet re-scored" : ""))));
  head.appendChild(ctlRow);

  // Action buttons
  const actions = el("div", {style: "display:flex;gap:.4em;flex-wrap:wrap;margin-top:.4em"});
  actions.appendChild(el("button",
    {class: "btn primary small", disabled: !editCount,
     onclick: () => doRescore().then(() => renderTweak())},
    "▶ Re-score"));
  actions.appendChild(el("button",
    {class: "btn small", disabled: !editCount,
     onclick: () => { if (confirm("Discard all edits?")) { resetEdits(); renderTweak(); } }},
    "Reset edits"));
  actions.appendChild(el("button",
    {class: "btn small", disabled: !editCount,
     onclick: () => downloadTweakFile(mode, "fasta")},
    "⤓ Tweaked FASTA"));
  actions.appendChild(el("button",
    {class: "btn small", disabled: !editCount,
     onclick: () => downloadTweakFile(mode, "xlsx")},
    "⤓ Tweaked XLSX"));
  head.appendChild(actions);
  root.appendChild(head);

  // ── Edits list ──────────────────────────────────────────────────────────
  if (editCount === 0) {
    root.appendChild(el("div", {class: "section-card"},
      el("p", {class: "muted"},
        "No edits applied yet. Switch to ", el("strong", null, "Sequence View"),
        ", click ", el("strong", null, "Enable editing"),
        " in the sidebar, then click any residue in the ", label(mode), " row to edit it.")));
    return;
  }

  const editsCard = el("div", {class: "section-card"});
  editsCard.appendChild(el("h2", null, "Edits applied"));
  const editsTbl = el("table", {class: "metrics fit"});
  editsTbl.appendChild(el("thead", null,
    el("tr", null,
      el("th", null, "Chain"),
      el("th", null, "Linear pos"),
      el("th", null, "IMGT pos"),
      el("th", null, "Original"),
      el("th", null, "Edited"),
    )));
  const editsBody = el("tbody");
  for (const chain of ["VH", "VL"]) {
    const ed = (EDITS[mode] || {})[chain] || {};
    const imap = imgtMapFor(mode, chain) || [];
    const orig = seqFor(mode, chain);
    const keys = Object.keys(ed).map(Number).sort((a,b) => a-b);
    for (const k of keys) {
      const i = k - 1;
      editsBody.appendChild(el("tr", null,
        el("td", null, chain),
        el("td", null, k),
        el("td", null, imap[i] != null ? imap[i] : "?"),
        el("td", null, el("span", {class: "res res-empty"}, orig[i] || "?")),
        el("td", null, el("span", {class: "res res-back_mutated"}, ed[k])),
      ));
    }
  }
  editsTbl.appendChild(editsBody);
  editsCard.appendChild(editsTbl);
  root.appendChild(editsCard);

  // ── Comparison table: metric / Original / Tweaked / Δ ──────────────────
  if (!tw) {
    root.appendChild(el("div", {class: "section-card"},
      el("p", {class: "muted"},
        "Click ", el("strong", null, "▶ Re-score"),
        " to recompute OASis, CamSol, germline identity and structure (if enabled at submit time).")));
    return;
  }
  const orig = D.modes[mode].scores;

  const cmpCard = el("div", {class: "section-card"});
  cmpCard.appendChild(el("h2", null, "Original vs Tweaked"));
  const cmpTbl = el("table", {class: "metrics"});
  cmpTbl.appendChild(el("thead", null,
    el("tr", null,
      el("th", null, "Metric"),
      el("th", {class: "num"}, "Original"),
      el("th", {class: "num"}, "Tweaked"),
      el("th", {class: "num"}, "Δ"),
    )));
  const cmpBody = el("tbody");

  function group(name) {
    cmpBody.appendChild(el("tr", {class: "group-row"},
      el("td", {colspan: 4}, "─── " + name + " ───")));
  }
  function row(label, key, fmt) {
    const o = key.includes(".") ? key.split(".").reduce((a,k)=>a&&a[k], orig) : orig[key];
    const t = key.includes(".") ? key.split(".").reduce((a,k)=>a&&a[k], tw) : tw[key];
    cmpBody.appendChild(el("tr", null,
      el("td", null, label),
      el("td", {class: "num"}, o == null ? "—" : fmt(o)),
      el("td", {class: "num"}, t == null ? "—" : fmt(t)),
      el("td", {class: "num"}, fmtDelta(o, t, fmt)),
    ));
  }

  group("Humanness — OASis");
  row("OASis FR VH",       "oa_fr_vh", pct);
  row("OASis CDR VH",      "oa_cdr_vh", pct);
  row("Germline FR id VH", "g_vh", pct);
  row("OASis FR VL",       "oa_fr_vl", pct);
  row("OASis CDR VL",      "oa_cdr_vl", pct);
  row("Germline FR id VL", "g_vl", pct);

  group("Solubility — CamSol");
  row("CamSol FR VH",   "cs_fr_vh", v => num(v, 3));
  row("CamSol CDR VH",  "cs_cdr_vh", v => num(v, 3));
  row("Hotspots VH",    "hs_vh", intval);
  row("Hotspots FR VH", "hs_fr_vh", intval);
  row("CamSol FR VL",   "cs_fr_vl", v => num(v, 3));
  row("CamSol CDR VL",  "cs_cdr_vl", v => num(v, 3));
  row("Hotspots VL",    "hs_vl", intval);
  row("Hotspots FR VL", "hs_fr_vl", intval);

  group("Vernier / back-mutations");
  row("Vernier mutable VH",    "vh_vern_mut", intval);
  row("Vernier back-mut VH",   "vh_vern_back", intval);
  row("Vernier humanized VH",  "vh_vern_hum", intval);
  row("Vernier mutable VL",    "vl_vern_mut", intval);
  row("Vernier back-mut VL",   "vl_vern_back", intval);
  row("Vernier humanized VL",  "vl_vern_hum", intval);

  group("FR liabilities");
  row("Deamid FR VH",  "lia.vh_df", intval);
  row("Oxidation FR VH","lia.vh_of", intval);
  row("Isomer FR VH",  "lia.vh_if", intval);
  row("Deamid FR VL",  "lia.vl_df", intval);
  row("Oxidation FR VL","lia.vl_of", intval);
  row("Isomer FR VL",  "lia.vl_if", intval);

  group("Physicochemical");
  row("pI Fv",         "pi", v => num(v, 2));
  row("Net charge Fv", "ch", v => num(v, 2));

  const anyStruct = tw.cf_vh != null || orig.cf_vh != null;
  if (anyStruct) {
    group("Structure (ABodyBuilder2)");
    row("Conf mean VH",  "cf_vh", v => num(v, 3));
    row("Conf FR VH",    "cf_fr_vh", v => num(v, 3));
    row("Conf CDR VH",   "cf_cdr_vh", v => num(v, 3));
    row("Conf CDR-H3",   "cf_cdr3_vh", v => num(v, 3));
    row("Conf min VH",   "cf_min_vh", v => num(v, 3));
    row("Conf mean VL",  "cf_vl", v => num(v, 3));
    row("Conf FR VL",    "cf_fr_vl", v => num(v, 3));
    row("Conf CDR VL",   "cf_cdr_vl", v => num(v, 3));
    row("Conf min VL",   "cf_min_vl", v => num(v, 3));
  }
  cmpTbl.appendChild(cmpBody);
  cmpCard.appendChild(cmpTbl);
  root.appendChild(cmpCard);
}

function downloadTweakFile(mode, format) {
  const body = {
    mode,
    edits: { VH: (EDITS[mode]||{}).VH || {}, VL: (EDITS[mode]||{}).VL || {} },
  };
  // Server applies edits on its side, re-scores and builds artifact
  fetch(`/api/report/${JOB}/tweak/${format}`, {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify(body),
  }).then(async resp => {
    if (!resp.ok) {
      alert("Download failed: " + (await resp.text()));
      return;
    }
    const blob = await resp.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `tweaked_${mode}_${JOB}.${format}`;
    a.click();
    URL.revokeObjectURL(url);
  });
}

// =====================================================================
// Tweak / rescore
// =====================================================================
function openEditPopup(ev, mode, chain, linIdx, imgtPos, currentAA) {
  // remove any existing popup
  document.querySelectorAll(".popup").forEach(p => p.remove());

  const status = (() => {
    const list = mutableFor(mode, chain);
    const hit = list.find(r => {
      // find by imgt
      return imgtPos != null && r.pos === imgtPos;
    });
    return hit ? hit.status : "—";
  })();

  const popup = el("div", {class: "popup"},
    el("div", null, el("strong", null, `IMGT pos ${imgtPos ?? "?"}`), " — current ", el("strong", null, currentAA),
      " (", status, ")"),
    el("label", null, "New residue:",
      el("input", {type: "text", id: "popup-input", maxlength: 1, value: currentAA, autocomplete: "off"})),
    el("div", {class: "actions"},
      el("button", {class: "btn small", onclick: () => popup.remove()}, "Cancel"),
      el("button", {class: "btn primary small", onclick: () => {
        const v = document.getElementById("popup-input").value.trim().toUpperCase();
        if (!v || !"ACDEFGHIKLMNPQRSTVWY".includes(v)) {
          alert("Invalid amino acid. Use one of: A C D E F G H I K L M N P Q R S T V W Y");
          return;
        }
        // record edit (or remove if same as original)
        const orig = D.modes[mode][chain === "VH" ? "vh" : "vl"][linIdx - 1];
        if (v === orig) {
          delete EDITS[mode][chain][linIdx];
        } else {
          EDITS[mode][chain][linIdx] = v;
        }
        popup.remove();
        updateRescoreBar();
        renderSequenceView();
      }}, "Apply"))
  );
  document.body.appendChild(popup);
  const rect = ev.target.getBoundingClientRect();
  popup.style.left = (rect.left + window.scrollX) + "px";
  popup.style.top  = (rect.bottom + window.scrollY + 4) + "px";
  setTimeout(() => document.getElementById("popup-input").focus(), 0);
}

function totalEdits() {
  let n = 0;
  for (const m of Object.keys(EDITS)) {
    for (const c of ["VH", "VL"]) n += Object.keys(EDITS[m][c] || {}).length;
  }
  return n;
}
function updateRescoreBar() {
  const bar = document.getElementById("rescore-bar");
  const n = totalEdits();
  updateTweakCountBadge();
  if (n === 0) {
    bar.hidden = true;
    return;
  }
  bar.hidden = false;
  // Summary: list edits per mode/chain
  const parts = [];
  for (const m of Object.keys(EDITS)) {
    for (const c of ["VH", "VL"]) {
      const ed = EDITS[m][c];
      const keys = Object.keys(ed);
      if (!keys.length) continue;
      const pieces = keys.map(k => {
        const orig = D.modes[m][c === "VH" ? "vh" : "vl"][+k - 1];
        return `${k}(${orig}→${ed[k]})`;
      });
      parts.push(`${label(m)}/${c}: ${pieces.join(",")}`);
    }
  }
  document.getElementById("rescore-summary").textContent = `${n} edit${n>1?"s":""}: ${parts.join(" · ")}`;
}

async function doRescore() {
  const goBtn = document.getElementById("rescore-go");
  goBtn.disabled = true;
  const labelOriginal = "▶ Re-score";
  try {
    // One request per mode, sending BOTH chains' edits at once. The server
    // also runs ABodyBuilder2 paired-Fv structure prediction if the original
    // run had structure enabled — slow (~30–60s) but acceptable for tweak.
    for (const m of Object.keys(EDITS)) {
      const vhEdits = EDITS[m].VH || {};
      const vlEdits = EDITS[m].VL || {};
      if (!Object.keys(vhEdits).length && !Object.keys(vlEdits).length) continue;
      goBtn.textContent = `Re-scoring ${label(m)}…`;
      const resp = await fetch(`/api/report/${JOB}/rescore`, {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({mode: m, edits: {VH: vhEdits, VL: vlEdits}}),
      });
      const j = await resp.json();
      if (j.error) { alert(`Rescore failed for ${label(m)}: ${j.error}`); continue; }
      if (!RESCORED[m]) RESCORED[m] = {};
      if (j.VH) RESCORED[m].VH = j.VH;
      if (j.VL) RESCORED[m].VL = j.VL;
      RESCORED[m].structure = j.structure || null;

      // Patch _all so scoresFor(m) returns merged values
      const base = Object.assign({}, D.modes[m].scores);
      base.lia = Object.assign({}, D.modes[m].scores.lia || {});
      for (const chain2 of ["VH", "VL"]) {
        const d2 = RESCORED[m][chain2];
        if (!d2) continue;
        if (chain2 === "VH") {
          base.oa_vh = d2.oasis; base.oa_fr_vh = d2.oasis_fr; base.oa_cdr_vh = d2.oasis_cdr;
          base.g_vh = d2.germline_id;
          base.cs_vh = d2.camsol; base.cs_fr_vh = d2.camsol_fr; base.cs_cdr_vh = d2.camsol_cdr;
          // Hotspot total = FR + CDR (server doesn't return the total separately)
          base.hs_vh = (d2.hs_fr ?? 0) + (d2.hs_cdr ?? 0);
          base.hs_fr_vh = d2.hs_fr; base.hs_cdr_vh = d2.hs_cdr;
          base.vh_vern_mut = d2.vern_mut; base.vh_vern_back = d2.vern_back;
          base.vh_vern_hum = d2.vern_hum; base.vh_vern_other = d2.vern_other;
          base.vh_fr_mut = d2.fr_mut; base.vh_fr_back = d2.fr_back;
          base.vh_fr_hum = d2.fr_hum; base.vh_fr_other = d2.fr_other;
          base.lia.vh_dc = d2.lia.dc; base.lia.vh_df = d2.lia.df;
          base.lia.vh_oc = d2.lia.oc; base.lia.vh_of = d2.lia.of;
          base.lia.vh_if = d2.lia.if;
        } else {
          base.oa_vl = d2.oasis; base.oa_fr_vl = d2.oasis_fr; base.oa_cdr_vl = d2.oasis_cdr;
          base.g_vl = d2.germline_id;
          base.cs_vl = d2.camsol; base.cs_fr_vl = d2.camsol_fr; base.cs_cdr_vl = d2.camsol_cdr;
          base.hs_vl = (d2.hs_fr ?? 0) + (d2.hs_cdr ?? 0);
          base.hs_fr_vl = d2.hs_fr; base.hs_cdr_vl = d2.hs_cdr;
          base.vl_vern_mut = d2.vern_mut; base.vl_vern_back = d2.vern_back;
          base.vl_vern_hum = d2.vern_hum; base.vl_vern_other = d2.vern_other;
          base.vl_fr_mut = d2.fr_mut; base.vl_fr_back = d2.fr_back;
          base.vl_fr_hum = d2.fr_hum; base.vl_fr_other = d2.fr_other;
          base.lia.vl_dc = d2.lia.dc; base.lia.vl_df = d2.lia.df;
          base.lia.vl_oc = d2.lia.oc; base.lia.vl_of = d2.lia.of;
          base.lia.vl_if = d2.lia.if;
        }
      }
      // Patch structure if the server returned new values
      const st = RESCORED[m].structure || {};
      const stKeys = [
        ["cf_vh","conf_mean_vh"], ["cf_vl","conf_mean_vl"],
        ["cf_fr_vh","conf_fr_mean_vh"], ["cf_fr_vl","conf_fr_mean_vl"],
        ["cf_cdr_vh","conf_cdr_mean_vh"], ["cf_cdr_vl","conf_cdr_mean_vl"],
        ["cf_cdr1_vh","conf_cdr1_mean_vh"], ["cf_cdr2_vh","conf_cdr2_mean_vh"],
        ["cf_cdr3_vh","conf_cdr3_mean_vh"],
        ["cf_cdr1_vl","conf_cdr1_mean_vl"], ["cf_cdr2_vl","conf_cdr2_mean_vl"],
        ["cf_cdr3_vl","conf_cdr3_mean_vl"],
        ["cf_min_vh","conf_min_vh"], ["cf_min_vl","conf_min_vl"],
      ];
      for (const [shortK, longK] of stKeys) {
        if (st[longK] != null) base[shortK] = +(+st[longK]).toFixed(4);
      }
      RESCORED[m]._all = base;
    }
    renderAll();
  } finally {
    goBtn.disabled = false;
    goBtn.textContent = labelOriginal;
  }
}

function resetEdits() {
  for (const m of Object.keys(EDITS)) {
    EDITS[m] = { VH: {}, VL: {} };
    delete RESCORED[m];
  }
  updateRescoreBar();
  renderAll();
}

function downloadTweakedFasta() {
  const lines = [];
  for (const m of sortedActive()) {
    if (m === "sapiens" || m === "lab") continue;
    for (const c of ["VH", "VL"]) {
      const orig = D.modes[m][c === "VH" ? "vh" : "vl"] || "";
      const ed = EDITS[m] && EDITS[m][c] || {};
      const chars = orig.split("");
      for (const [k, v] of Object.entries(ed)) {
        const i = +k - 1;
        if (i >= 0 && i < chars.length) chars[i] = v;
      }
      lines.push(`>${m}_${c}_${JOB}_tweaked`);
      lines.push(chars.join(""));
    }
  }
  const blob = new Blob([lines.join("\n") + "\n"], {type: "text/plain"});
  const a = el("a", {href: URL.createObjectURL(blob), download: `tweaked_${JOB}.fasta`});
  a.click();
}

// =====================================================================
// Init
// =====================================================================
function renderAll() {
  renderLegend();
  renderSummary();
  renderSequenceView();
  renderGermline();
  renderMutable();
  renderFeatures();
  renderTweak();
  updateTweakCountBadge();
}

function updateTweakCountBadge() {
  const badge = document.getElementById("tweak-count");
  if (!badge) return;
  const n = totalEdits();
  badge.textContent = n;
  badge.hidden = n === 0;
}

function initTabs() {
  document.querySelectorAll(".nav-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      const target = btn.dataset.tab;
      document.querySelectorAll(".nav-btn").forEach(b => b.classList.toggle("active", b === btn));
      document.querySelectorAll(".tab").forEach(t => t.hidden = t.dataset.tab !== target);
    });
  });
}

function initTweak() {
  const toggle = document.getElementById("tweak-toggle");
  const controls = document.getElementById("tweak-controls");
  const select = document.getElementById("tweak-mode");
  select.innerHTML = "";
  for (const m of sortedActive()) {
    if (m === "sapiens" || m === "lab") continue;  // not editable per spec §10.5
    const opt = el("option", {value: m}, label(m));
    select.appendChild(opt);
  }
  toggle.addEventListener("click", () => {
    window._TWEAK_ON = !window._TWEAK_ON;
    toggle.textContent = window._TWEAK_ON ? "Disable editing" : "Enable editing";
    controls.hidden = !window._TWEAK_ON;
    renderSequenceView();
  });

  document.getElementById("rescore-go").addEventListener("click", doRescore);
  document.getElementById("rescore-reset").addEventListener("click", () => {
    if (confirm("Discard all edits?")) resetEdits();
  });
  document.getElementById("rescore-download").addEventListener("click", downloadTweakedFasta);
}

function initDownloads() {
  document.getElementById("dl-xlsx").addEventListener("click", () => {
    window.location = `/api/report/${JOB}/xlsx`;
  });
  document.getElementById("dl-html").addEventListener("click", () => {
    window.location = `/api/report/${JOB}/html`;
  });
  document.getElementById("dl-fasta").addEventListener("click", () => {
    window.location = `/api/report/${JOB}/fasta`;
  });
}

// Close any open popup on click outside
document.addEventListener("click", (ev) => {
  if (ev.target.closest(".popup") || ev.target.closest(".aa.tweak-clickable")) return;
  document.querySelectorAll(".popup").forEach(p => p.remove());
});

initTabs();
initTweak();
initDownloads();
renderAll();
})();
