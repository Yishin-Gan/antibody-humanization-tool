// Antibody Humanization Advisor — input page

const VALID_AA = new Set("ACDEFGHIKLMNPQRSTVWY");

function cleanAA(raw) {
  if (!raw) return "";
  const out = [];
  for (const line of raw.split(/\r?\n/)) {
    if (line.startsWith(">")) continue;
    for (const ch of line) {
      if (/[a-z]/i.test(ch)) {
        const u = ch.toUpperCase();
        if (VALID_AA.has(u)) out.push(u);
      }
    }
  }
  return out.join("");
}

function setHint(el, msg, cls) {
  el.textContent = msg || "";
  el.className = "hint" + (cls ? " " + cls : "");
}

// VH/VL paste hint
function aaCount(textarea, hint, max, label) {
  textarea.addEventListener("input", () => {
    const cleaned = cleanAA(textarea.value);
    if (!cleaned.length) {
      setHint(hint, "");
      return;
    }
    if (cleaned.length < 80) {
      setHint(hint, `${cleaned.length} aa — too short (min 80)`, "error");
    } else if (cleaned.length > max) {
      setHint(hint, `${cleaned.length} aa — too long (max ${max})`, "warn");
    } else {
      setHint(hint, `${cleaned.length} aa — OK`, "ok");
    }
  });
}

aaCount(document.getElementById("mouse_vh"), document.getElementById("vh-hint"), 150, "VH");
aaCount(document.getElementById("mouse_vl"), document.getElementById("vl-hint"), 130, "VL");

// VL chain-type auto-detect
const vlField = document.getElementById("mouse_vl");
const detectEl = document.getElementById("chain-detect");
let detectTimer = null;
let detectedChain = null;
vlField.addEventListener("blur", async () => {
  const seq = cleanAA(vlField.value);
  if (seq.length < 60) return;
  setHint(detectEl, "detecting…");
  try {
    const r = await fetch("/api/detect-vl-chain", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({seq})
    });
    const j = await r.json();
    if (j.ok && j.chain_type) {
      detectedChain = j.chain_type;
      setHint(detectEl, `Auto-detected: ${j.chain_type === "K" ? "Kappa" : "Lambda"} ✓`, "ok");
      // If user has Auto selected, set radio to detected chain when submitting
    } else {
      setHint(detectEl, "Auto-detection uncertain — please select", "warn");
    }
  } catch (e) {
    setHint(detectEl, "detection failed: " + e, "error");
  }
});

// Preferred germline validation
async function validateGermline(input, hintEl, chainHint) {
  const name = input.value.trim();
  if (!name) { setHint(hintEl, ""); return; }
  let chain = chainHint;
  if (chain === "auto") {
    chain = detectedChain || (document.querySelector('input[name=vl_chain_type]:checked').value !== "auto"
                              ? document.querySelector('input[name=vl_chain_type]:checked').value
                              : "K");
  }
  try {
    const r = await fetch("/api/validate-germline", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({name, chain})
    });
    const j = await r.json();
    if (j.found) {
      setHint(hintEl, `✓ ${j.resolved}${j.note ? " (" + j.note + ")" : ""}`, "ok");
    } else {
      setHint(hintEl, `⚠ Not found in database — germline identity and FR alignment will be unavailable`, "warn");
    }
  } catch (e) {
    setHint(hintEl, "validation failed: " + e, "error");
  }
}

document.querySelectorAll("[data-validate-germline]").forEach(input => {
  const hintId = input.id.replace(/_/g, "-").replace(/germ/, "").replace(/^pref/, "pref")  + "-hint";
  // Map ids: preferred_germ_vh -> pref-vh-hint, lab_germ_vh -> lab-vh-hint
  let hintEl;
  if (input.id.startsWith("preferred_germ_")) {
    hintEl = document.getElementById("pref-" + input.id.split("_").pop() + "-hint");
  } else if (input.id.startsWith("lab_germ_")) {
    hintEl = document.getElementById("lab-" + input.id.split("_").pop() + "-hint");
  }
  if (!hintEl) return;
  const chain = input.dataset.validateGermline;
  input.addEventListener("blur", () => validateGermline(input, hintEl, chain));
});

// Lab section: auto-enable lab mode + require all-or-none
const labFields = [
  "lab_hu_vh", "lab_hu_vl",
  "lab_final_vh", "lab_final_vl",
  "lab_germ_vh", "lab_germ_vl",
];
const labModeBox = document.getElementById("mode-lab");
function updateLabMode() {
  const anyFilled = labFields.some(id => document.getElementById(id).value.trim().length > 0);
  labModeBox.checked = anyFilled;
}
labFields.forEach(id => {
  document.getElementById(id).addEventListener("input", updateLabMode);
});

// Submit handler: validation + progress overlay
const form = document.getElementById("run-form");
const runBtn = document.getElementById("run-btn");
const progress = document.getElementById("progress");
form.addEventListener("submit", async (e) => {
  e.preventDefault();  // we POST via fetch so the overlay stays visible
  // basic validation
  const vh = cleanAA(document.getElementById("mouse_vh").value);
  const vl = cleanAA(document.getElementById("mouse_vl").value);
  if (vh.length < 80 || vh.length > 150) {
    alert(`Mouse VH must be 80–150 amino acids (got ${vh.length}).`);
    return;
  }
  if (vl.length < 80 || vl.length > 130) {
    alert(`Mouse VL must be 80–130 amino acids (got ${vl.length}).`);
    return;
  }
  const selectedModes = Array.from(document.querySelectorAll('input[name="modes"]:checked')).map(c => c.value);
  if (selectedModes.length === 0) {
    alert("Select at least one processing mode.");
    return;
  }
  // Lab-mode all-or-none enforcement
  if (labModeBox.checked) {
    const missing = labFields.filter(id => !document.getElementById(id).value.trim());
    if (missing.length) {
      alert("Lab reference: all six fields required.\nMissing: " + missing.join(", "));
      return;
    }
  }
  // Preferred mode: at least one germline supplied
  const prefBox = document.getElementById("mode-preferred");
  if (prefBox.checked) {
    const pvh = document.getElementById("preferred_germ_vh").value.trim();
    const pvl = document.getElementById("preferred_germ_vl").value.trim();
    if (!pvh && !pvl) {
      alert("Preferred germline mode: supply at least one of VH/VL germline names.");
      return;
    }
  }

  // Build JSON body — preserves array semantics for `modes`.
  const labCleaned = {};
  for (const id of ["lab_hu_vh", "lab_hu_vl", "lab_final_vh", "lab_final_vl"]) {
    labCleaned[id] = cleanAA(document.getElementById(id).value);
  }
  const body = {
    mouse_vh: vh,
    mouse_vl: vl,
    vl_chain_type: document.querySelector('input[name="vl_chain_type"]:checked').value,
    modes: selectedModes,
    preferred_germ_vh: document.getElementById("preferred_germ_vh").value.trim(),
    preferred_germ_vl: document.getElementById("preferred_germ_vl").value.trim(),
    lab_hu_vh: labCleaned.lab_hu_vh,
    lab_hu_vl: labCleaned.lab_hu_vl,
    lab_final_vh: labCleaned.lab_final_vh,
    lab_final_vl: labCleaned.lab_final_vl,
    lab_germ_vh: document.getElementById("lab_germ_vh").value.trim(),
    lab_germ_vl: document.getElementById("lab_germ_vl").value.trim(),
    structure: document.getElementById("structure").checked ? "on" : "",
  };

  runBtn.disabled = true;
  runBtn.textContent = "Running…";
  progress.hidden = false;

  // Fake stage progression so user sees motion (real run takes 30s-2min)
  const stages = ["numbering", "germline", "graft", "sapiens", "score", "structure"];
  let i = 0;
  let tickHandle = null;
  function tick() {
    if (i > 0) {
      const prev = document.querySelector(`.stage[data-stage="${stages[i-1]}"]`);
      if (prev) { prev.classList.remove("active"); prev.classList.add("done");
                  prev.textContent = "✓ " + prev.textContent.replace(/^[…⏳✓]\s*/, ""); }
    }
    if (i < stages.length) {
      const cur = document.querySelector(`.stage[data-stage="${stages[i]}"]`);
      if (cur) { cur.classList.add("active");
                 cur.textContent = "⏳ " + cur.textContent.replace(/^[…⏳✓]\s*/, ""); }
      i++;
      tickHandle = setTimeout(tick, 6000);
    }
  }
  tick();

  console.log("[submit] starting fetch /run", body);
  try {
    const resp = await fetch("/run", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify(body),
    });
    if (tickHandle) clearTimeout(tickHandle);
    console.log("[submit] response:", resp.status, "redirected:", resp.redirected, "url:", resp.url,
                "content-type:", resp.headers.get("content-type"));

    document.querySelectorAll(".stage").forEach(s => {
      s.classList.remove("active"); s.classList.add("done");
      s.textContent = "✓ " + s.textContent.replace(/^[…⏳✓]\s*/, "");
    });

    // Case 1: server returned a 302 redirect → fetch followed it → final URL = report URL
    if (resp.redirected && resp.url && /\/report\//.test(resp.url)) {
      console.log("[submit] navigating to (redirected):", resp.url);
      window.location.assign(resp.url);
      return;
    }

    // Case 2: server returned JSON with {job_id, report_url} or {error}
    let j = null;
    const ct = (resp.headers.get("content-type") || "").toLowerCase();
    if (ct.includes("application/json")) {
      j = await resp.json();
    } else {
      const txt = await resp.text();
      try { j = JSON.parse(txt); }
      catch { throw new Error(txt.slice(0, 300) || `Server returned ${resp.status}`); }
    }
    console.log("[submit] parsed json:", j);
    if (!resp.ok || j.error) throw new Error(j.error || `Server returned ${resp.status}`);
    const target = j.report_url || (j.job_id ? "/report/" + j.job_id : null);
    if (!target) throw new Error("Server returned no report URL: " + JSON.stringify(j));
    console.log("[submit] navigating to:", target);

    // Show a fallback link inside the overlay so the user can recover if
    // the browser blocks programmatic navigation for any reason.
    const box = document.querySelector(".progress-box");
    if (box) {
      const link = document.createElement("p");
      link.style.marginTop = "1em";
      link.innerHTML = `Done. <a href="${target}" style="color:var(--accent)">Open report →</a>`;
      box.appendChild(link);
    }
    // Programmatic navigation
    window.location.assign(target);
    // Belt-and-braces: in case assign() is silently no-op, also try setting href after a tick
    setTimeout(() => { if (window.location.pathname === "/") window.location.href = target; }, 200);
  } catch (err) {
    if (tickHandle) clearTimeout(tickHandle);
    console.error("[submit] failed:", err);
    progress.hidden = true;
    runBtn.disabled = false;
    runBtn.textContent = "▶ Run humanization";
    alert("Run failed: " + err.message);
  }
});
