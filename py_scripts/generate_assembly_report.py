#!/usr/bin/env python3
"""
generate_assembly_report.py
===========================
Parse assembly QC metrics CSV and generate a self-contained interactive HTML
report that ranks assemblies by annotation suitability.

Ranking criteria (weighted, user-adjustable in report):
  - Contiguity:   N50 (scaffold-level preferred), L50, # contigs
  - Base accuracy: Merqury QV
  - Gene space:   BUSCO complete %
  - Scaffolding:  Hi-C trans/cis ratio (lower = better)
  - Completeness: Diploid k-mer completeness

Usage:
    python3 generate_assembly_report.py \
        --input assembly_qc_metrics.csv \
        --output assembly_qc_report.html \
        --report-stage final
"""

import argparse
import csv
import json
import sys
from collections import defaultdict

# ── Metric extraction helpers ────────────────────────────────────────────────

METRIC_MAP = {
    "quast": {
        "N50":                      {"key": "n50",      "scale": 1e-6, "unit": "Mb"},
        "L50":                      {"key": "l50",      "scale": 1,    "unit": ""},
        "# contigs":                {"key": "contigs",  "scale": 1,    "unit": ""},
        "Largest contig":           {"key": "largest",  "scale": 1e-6, "unit": "Mb"},
        "Total length":             {"key": "total",    "scale": 1e-6, "unit": "Mb"},
        "GC (%)":                   {"key": "gc",       "scale": 1,    "unit": "%"},
        "auN":                      {"key": "aun",      "scale": 1e-6, "unit": "Mb"},
        "N90":                      {"key": "n90",      "scale": 1e-6, "unit": "Mb"},
        "# N's per 100 kbp":       {"key": "ns_100k",  "scale": 1,    "unit": ""},
    },
    "busco": {
        "complete":     {"key": "busco_c",    "scale": 1, "unit": ""},
        "single":       {"key": "busco_s",    "scale": 1, "unit": ""},
        "duplicated":   {"key": "busco_d",    "scale": 1, "unit": ""},
        "fragmented":   {"key": "busco_f",    "scale": 1, "unit": ""},
        "missing":      {"key": "busco_m",    "scale": 1, "unit": ""},
        "total_busco":  {"key": "busco_total","scale": 1, "unit": ""},
    },
    "merqury": {
        "qv":                {"key": "qv",          "scale": 1, "unit": ""},
        "kmer_completeness": {"key": "kcomp",       "scale": 1, "unit": "%"},
        "error_kmer":        {"key": "error_kmer",  "scale": 1, "unit": ""},
    },
    "hic_contact": {
        "trans_to_cis":   {"key": "tc_ratio",      "scale": 1, "unit": ""},
        "cis_pairs":      {"key": "cis_pairs",     "scale": 1, "unit": ""},
        "trans_pairs":    {"key": "trans_pairs",    "scale": 1, "unit": ""},
        "retention_pct":  {"key": "retention_pct",  "scale": 1, "unit": "%"},
    },
    "mapped_hifi": {
        "hifi_depth":  {"key": "hifi_depth", "scale": 1, "unit": "x"},
    },
}


def parse_metrics(csv_path):
    """Parse the long-format metrics CSV into a nested dict structure."""
    data = defaultdict(lambda: defaultdict(lambda: {"hap1": {}, "hap2": {}, "both": {}}))

    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid   = row["sample_id"]
            stage = row["stage"]
            analysis = row["analysis"]
            metric   = row["metric"]

            if analysis not in METRIC_MAP or metric not in METRIC_MAP[analysis]:
                continue

            info = METRIC_MAP[analysis][metric]
            key  = info["key"]
            scale = info["scale"]

            for hap in ["hap1", "hap2", "both"]:
                raw = row.get(hap, "").strip()
                if raw and raw != "NA":
                    try:
                        data[sid][stage][hap][key] = float(raw) * scale
                    except ValueError:
                        pass

    return data


def get_available_stages(data):
    """Return ordered list of all stages present in the data."""
    stages = set()
    for sid in data:
        stages.update(data[sid].keys())

    # Full pipeline order — must match labels from compile_qc.R
    stage_order = [
        "ctg.base", "ctg.mito", "ctg.purged", "ctg.cor", "ctg.deco",
        "scaf.base", "scaf.cor", "scaf2",
        "gap_fill", "final"
    ]
    ordered = [s for s in stage_order if s in stages]
    remaining = sorted(stages - set(ordered))
    return ordered + remaining


def build_report_data(data, stages):
    """Build the JSON data structure for the HTML report."""
    report = {
        "stages": stages,
        "samples": sorted(data.keys()),
        "stage_data": {},
    }

    for stage in stages:
        stage_info = {"assemblies": {}}

        for sid in data:
            if stage not in data[sid]:
                continue

            rec = data[sid][stage]
            assembly = {}
            for hap in ["hap1", "hap2", "both"]:
                assembly[hap] = {}
                for k, v in rec[hap].items():
                    assembly[hap][k] = round(v, 4) if isinstance(v, float) else v

            # Add computed fields
            for hap in ["hap1", "hap2"]:
                h = assembly[hap]
                if "busco_c" in h and "busco_total" in h and h["busco_total"] > 0:
                    h["busco_pct"] = round(h["busco_c"] / h["busco_total"] * 100, 2)

            stage_info["assemblies"][sid] = assembly

        report["stage_data"][stage] = stage_info

    return report


# ── HTML template ────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Assembly QC Report</title>
<style>
  :root {
    --bg: #0f172a; --surface: #1e293b; --border: #334155;
    --text: #e2e8f0; --muted: #94a3b8; --dim: #64748b;
    --green: #22c55e; --blue: #3b82f6; --purple: #a855f7;
    --amber: #f59e0b; --red: #ef4444; --cyan: #06b6d4;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: var(--bg); color: var(--text); font-family: system-ui, -apple-system, sans-serif; padding: 20px; }
  .container { max-width: 1100px; margin: 0 auto; }
  h1 { font-size: 24px; font-weight: 700; margin-bottom: 4px; }
  .subtitle { font-size: 13px; color: var(--dim); margin-bottom: 20px; }

  /* Controls */
  .controls { display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 16px; align-items: center; }
  .controls label { font-size: 12px; color: var(--muted); font-weight: 600; margin-right: 4px; }
  .btn { padding: 7px 14px; border-radius: 8px; border: 2px solid var(--border); background: var(--surface);
         color: var(--muted); cursor: pointer; font-size: 13px; font-weight: 600; transition: all 0.2s; }
  .btn:hover { border-color: var(--blue); color: var(--blue); }
  .btn.active { border-color: var(--blue); background: rgba(59,130,246,0.12); color: var(--blue); }

  /* Weights panel */
  .weights-panel { background: var(--surface); border-radius: 12px; padding: 16px 20px; margin-bottom: 16px;
                   border: 1px solid var(--border); }
  .weights-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px; cursor: pointer; }
  .weights-header h3 { font-size: 13px; font-weight: 700; color: var(--muted); text-transform: uppercase; letter-spacing: 0.5px; }
  .weights-toggle { font-size: 11px; color: var(--dim); background: var(--bg); padding: 3px 10px; border-radius: 6px;
                    border: 1px solid var(--border); cursor: pointer; font-weight: 600; }
  .weights-toggle:hover { color: var(--blue); border-color: var(--blue); }
  .weights-body { display: none; }
  .weights-body.open { display: block; }
  .weights-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 12px 20px; }
  .weight-item { }
  .weight-label { display: flex; justify-content: space-between; align-items: center; margin-bottom: 4px; }
  .weight-name { font-size: 12px; color: var(--muted); font-weight: 600; }
  .weight-val { font-size: 12px; color: var(--blue); font-weight: 700; min-width: 28px; text-align: right; }
  .weight-slider { -webkit-appearance: none; appearance: none; width: 100%; height: 5px; border-radius: 3px;
                   background: var(--bg); outline: none; }
  .weight-slider::-webkit-slider-thumb { -webkit-appearance: none; appearance: none; width: 14px; height: 14px;
                                          border-radius: 50%; background: var(--blue); cursor: pointer; border: 2px solid var(--surface); }
  .weight-slider::-moz-range-thumb { width: 14px; height: 14px; border-radius: 50%; background: var(--blue);
                                      cursor: pointer; border: 2px solid var(--surface); }
  .weight-total { margin-top: 12px; display: flex; justify-content: space-between; align-items: center; }
  .weight-total-label { font-size: 12px; color: var(--dim); }
  .weight-total-val { font-size: 13px; font-weight: 700; }
  .weight-total-val.ok { color: var(--green); }
  .weight-total-val.warn { color: var(--amber); }
  .weight-actions { display: flex; gap: 8px; margin-top: 10px; }
  .weight-reset-btn { font-size: 11px; padding: 4px 12px; border-radius: 6px; border: 1px solid var(--border);
                      background: transparent; color: var(--muted); cursor: pointer; font-weight: 600; }
  .weight-reset-btn:hover { border-color: var(--amber); color: var(--amber); }
  .weight-preset-btn { font-size: 11px; padding: 4px 12px; border-radius: 6px; border: 1px solid var(--border);
                       background: transparent; color: var(--muted); cursor: pointer; font-weight: 600; }
  .weight-preset-btn:hover { border-color: var(--cyan); color: var(--cyan); }
  .weight-preset-btn.active-preset { border-color: var(--cyan); color: var(--cyan); background: rgba(6,182,212,0.1); }

  /* Cards and tables */
  .card { background: var(--surface); border-radius: 12px; padding: 20px; margin-bottom: 16px; border: 1px solid var(--border); }
  .card.selected { border-color: var(--sel-color, var(--blue)); }
  .badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px; }
  .sample-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px; }
  .sample-name { font-size: 20px; font-weight: 700; margin-top: 6px; }
  .metrics-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
  .metric-section-title { font-size: 11px; color: var(--dim); font-weight: 600; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 10px; }
  .metric-row { margin-bottom: 6px; }
  .metric-label { display: flex; justify-content: space-between; font-size: 12px; color: var(--muted); margin-bottom: 2px; }
  .metric-value { font-weight: 600; color: var(--text); }
  .bar-bg { height: 6px; background: var(--bg); border-radius: 3px; }
  .bar-fill { height: 6px; border-radius: 3px; transition: width 0.4s ease; }
  .summary-bar { margin-top: 14px; padding: 10px 14px; background: var(--bg); border-radius: 8px; font-size: 12px; color: var(--muted); }
  .summary-bar span.hl { font-weight: 600; }
  table { width: 100%; border-collapse: collapse; font-size: 12px; }
  th { padding: 6px 8px; text-align: left; color: var(--dim); font-weight: 600; border-bottom: 1px solid var(--border); }
  td { padding: 6px 8px; border-bottom: 1px solid rgba(51,65,85,0.4); }
  tr.clickable { cursor: pointer; transition: background 0.15s; }
  tr.clickable:hover { background: rgba(59,130,246,0.06); }
  tr.clickable.active { background: rgba(59,130,246,0.1); }
  .rank-cell { font-weight: 700; }
  .score-cell { font-weight: 600; }
  .tc-good { color: var(--green); } .tc-mid { color: var(--amber); } .tc-bad { color: var(--red); }
  .hap-btns { display: flex; gap: 4px; }
  .hap-btn { padding: 5px 12px; border-radius: 6px; border: 1px solid var(--border); background: transparent;
             color: var(--muted); cursor: pointer; font-size: 12px; font-weight: 600; transition: all 0.15s; }
  .hap-btn.active { border-color: var(--sel-color, var(--blue)); background: rgba(59,130,246,0.12); color: var(--sel-color, var(--blue)); }
  @media (max-width: 700px) { .metrics-grid { grid-template-columns: 1fr; } .weights-grid { grid-template-columns: 1fr 1fr; } }
</style>
</head>
<body>
<div class="container">
  <h1>Assembly QC Report</h1>
  <p class="subtitle">Interactive comparison of genome assemblies — ranked by annotation suitability</p>

  <div class="controls">
    <label>Stage:</label>
    <div id="stage-btns"></div>
  </div>

  <div class="weights-panel">
    <div class="weights-header" onclick="toggleWeights()">
      <h3>Ranking Weights</h3>
      <button class="weights-toggle" id="weights-toggle-btn">Show</button>
    </div>
    <div class="weights-body" id="weights-body">
      <div id="presets-row" class="weight-actions" style="margin-bottom:12px;margin-top:0"></div>
      <div class="weights-grid" id="weights-grid"></div>
      <div class="weight-total">
        <span class="weight-total-label">Total weight:</span>
        <span class="weight-total-val" id="weight-total">100</span>
      </div>
      <div class="weight-actions">
        <button class="weight-reset-btn" onclick="resetWeights()">Reset to defaults</button>
      </div>
    </div>
  </div>

  <div id="comparison-table" class="card"></div>
  <div id="detail-card"></div>
</div>

<script>
// ── Data injected by Python ──
const REPORT = __REPORT_JSON__;

const RANK_COLORS = ["#22c55e","#3b82f6","#a855f7","#f59e0b","#ef4444","#06b6d4","#fb7185","#a3e635","#818cf8","#fbbf24"];
const RANK_LABELS = ["BEST","STRONG","GOOD","MODERATE","BELOW AVG","POOR","POOR","POOR","POOR","POOR"];

// ── Weight definitions ──
const WEIGHT_DEFS = [
  { id: "n50",     name: "N50 (contiguity)",    default: 25, min: 0, max: 50 },
  { id: "qv",      name: "QV (base accuracy)",  default: 25, min: 0, max: 50 },
  { id: "busco",   name: "BUSCO complete",       default: 20, min: 0, max: 50 },
  { id: "tc",      name: "Trans/Cis ratio",      default: 15, min: 0, max: 50 },
  { id: "contigs", name: "Contig count",          default: 10, min: 0, max: 50 },
  { id: "kcomp",   name: "K-mer completeness",   default: 5,  min: 0, max: 50 },
];

const PRESETS = {
  "Default": { n50: 25, qv: 25, busco: 20, tc: 15, contigs: 10, kcomp: 5 },
  "Accuracy-first": { n50: 10, qv: 40, busco: 25, tc: 5, contigs: 5, kcomp: 15 },
  "Contiguity-first": { n50: 40, qv: 15, busco: 15, tc: 15, contigs: 10, kcomp: 5 },
  "Gene-space-first": { n50: 15, qv: 20, busco: 35, tc: 10, contigs: 10, kcomp: 10 },
  "Scaffolding-first": { n50: 15, qv: 15, busco: 15, tc: 35, contigs: 10, kcomp: 10 },
};

let weights = {};
WEIGHT_DEFS.forEach(w => weights[w.id] = w.default);

let curStage = REPORT.stages.includes("__DEFAULT_STAGE__") ? "__DEFAULT_STAGE__" : REPORT.stages[REPORT.stages.length - 1];
let curSample = null;
let curHap = "hap2";
let weightsOpen = false;

function $(id) { return document.getElementById(id); }

// ── Ranking engine (runs in browser with current weights) ──
function scoreHaplotype(h, w) {
  let s = 0, wt = 0;

  if (h.n50 !== undefined && w.n50 > 0) {
    const v = Math.max(0, Math.min(100, (Math.log10(Math.max(h.n50, 0.1)) - 0) / 2 * 100));
    s += v * w.n50; wt += w.n50;
  }
  if (h.qv !== undefined && w.qv > 0) {
    const v = Math.max(0, Math.min(100, (h.qv - 55) / 10 * 100));
    s += v * w.qv; wt += w.qv;
  }
  if (h.busco_c !== undefined && h.busco_total > 0 && w.busco > 0) {
    const pct = h.busco_c / h.busco_total * 100;
    const v = Math.max(0, Math.min(100, (pct - 85) / 10 * 100));
    s += v * w.busco; wt += w.busco;
  }
  if (h.contigs !== undefined && w.contigs > 0) {
    const v = Math.max(0, Math.min(100,
      (Math.log10(3000) - Math.log10(Math.max(h.contigs, 100))) / (Math.log10(3000) - Math.log10(100)) * 100));
    s += v * w.contigs; wt += w.contigs;
  }
  if (h.tc_ratio !== undefined && w.tc > 0) {
    const v = Math.max(0, Math.min(100, (1.8 - h.tc_ratio) / 1.5 * 100));
    s += v * w.tc; wt += w.tc;
  }
  if (h.kcomp !== undefined && w.kcomp > 0) {
    const v = Math.max(0, Math.min(100, (h.kcomp - 70) / 12 * 100));
    s += v * w.kcomp; wt += w.kcomp;
  }

  return wt > 0 ? s / wt : 0;
}

function rankAssemblies(sd) {
  const results = [];
  for (const sid of Object.keys(sd.assemblies)) {
    const a = sd.assemblies[sid];
    const hapScores = {};
    for (const hap of ["hap1","hap2"]) {
      const h = a[hap];
      if (h && Object.keys(h).length > 0) {
        hapScores[hap] = scoreHaplotype(h, weights);
      }
    }
    if (Object.keys(hapScores).length === 0) continue;
    const avg = Object.values(hapScores).reduce((a,b)=>a+b,0) / Object.values(hapScores).length;
    const bestHap = Object.entries(hapScores).sort((a,b)=>b[1]-a[1])[0][0];
    results.push({ sample_id: sid, score: avg, best_hap: bestHap });
  }
  results.sort((a,b) => b.score - a.score);
  results.forEach((r,i) => r.rank = i + 1);
  return results;
}

// ── Weights UI ──
function toggleWeights() {
  weightsOpen = !weightsOpen;
  $("weights-body").classList.toggle("open", weightsOpen);
  $("weights-toggle-btn").textContent = weightsOpen ? "Hide" : "Show";
}

function renderWeightsUI() {
  // Presets
  const pr = $("presets-row");
  pr.innerHTML = '<span style="font-size:11px;color:var(--dim);font-weight:600;margin-right:4px">Presets:</span>';
  for (const [name, vals] of Object.entries(PRESETS)) {
    const b = document.createElement("button");
    b.className = "weight-preset-btn";
    const isActive = WEIGHT_DEFS.every(w => weights[w.id] === vals[w.id]);
    if (isActive) b.classList.add("active-preset");
    b.textContent = name;
    b.onclick = () => { applyPreset(vals); };
    pr.appendChild(b);
  }

  // Sliders
  const grid = $("weights-grid");
  grid.innerHTML = "";
  WEIGHT_DEFS.forEach(wd => {
    const item = document.createElement("div");
    item.className = "weight-item";
    item.innerHTML = `
      <div class="weight-label">
        <span class="weight-name">${wd.name}</span>
        <span class="weight-val" id="wv-${wd.id}">${weights[wd.id]}</span>
      </div>
      <input type="range" class="weight-slider" id="ws-${wd.id}"
             min="${wd.min}" max="${wd.max}" value="${weights[wd.id]}" step="1">
    `;
    grid.appendChild(item);

    item.querySelector("input").addEventListener("input", function() {
      weights[wd.id] = parseInt(this.value);
      $("wv-" + wd.id).textContent = this.value;
      updateWeightTotal();
      renderAll();
      renderWeightsUI(); // refresh preset highlights
    });
  });

  updateWeightTotal();
}

function updateWeightTotal() {
  const total = WEIGHT_DEFS.reduce((s,w) => s + weights[w.id], 0);
  const el = $("weight-total");
  el.textContent = total;
  el.className = "weight-total-val " + (total === 100 ? "ok" : "warn");
}

function resetWeights() {
  applyPreset(PRESETS["Default"]);
}

function applyPreset(vals) {
  WEIGHT_DEFS.forEach(w => {
    weights[w.id] = vals[w.id];
    const slider = $("ws-" + w.id);
    if (slider) slider.value = vals[w.id];
    const display = $("wv-" + w.id);
    if (display) display.textContent = vals[w.id];
  });
  updateWeightTotal();
  renderAll();
  renderWeightsUI();
}

// ── Init ──
function init() {
  renderStageBtns();
  renderWeightsUI();
  renderAll();
}

function renderStageBtns() {
  const c = $("stage-btns");
  c.innerHTML = "";
  REPORT.stages.forEach(s => {
    const b = document.createElement("button");
    b.className = "btn" + (s === curStage ? " active" : "");
    b.textContent = s;
    b.onclick = () => { curStage = s; curSample = null; renderAll(); };
    c.appendChild(b);
  });
}

function renderAll() {
  document.querySelectorAll("#stage-btns .btn").forEach(b => {
    b.classList.toggle("active", b.textContent === curStage);
  });

  const sd = REPORT.stage_data[curStage];
  if (!sd) { $("comparison-table").innerHTML = "<p>No data for this stage.</p>"; $("detail-card").innerHTML = ""; return; }

  const ranking = rankAssemblies(sd);
  if (!curSample && ranking.length > 0) curSample = ranking[0].sample_id;

  renderTable(sd, ranking);
  renderDetail(sd, ranking);
}

function renderTable(sd, ranking) {
  const hasTc = ranking.some(r => {
    const a = sd.assemblies[r.sample_id];
    return a && (a.hap1.tc_ratio !== undefined || a.hap2.tc_ratio !== undefined);
  });

  let html = `<div class="metric-section-title">Assembly Ranking — ${curStage}</div>`;
  html += `<table><thead><tr><th>Rank</th><th>Sample</th><th>Score</th><th>N50 (Mb)</th><th>Contigs</th><th>BUSCO %</th><th>QV</th>`;
  if (hasTc) html += `<th>T/C Ratio</th>`;
  html += `<th>Best Hap</th></tr></thead><tbody>`;

  ranking.forEach((r, i) => {
    const a = sd.assemblies[r.sample_id];
    const color = RANK_COLORS[Math.min(i, RANK_COLORS.length - 1)];
    const h = a[r.best_hap] || {};
    const active = r.sample_id === curSample ? " active" : "";

    const tcVal = h.tc_ratio;
    const tcClass = tcVal !== undefined ? (tcVal < 0.5 ? "tc-good" : tcVal < 1 ? "tc-mid" : "tc-bad") : "";

    html += `<tr class="clickable${active}" onclick="selectSample('${r.sample_id}')">`;
    html += `<td class="rank-cell" style="color:${color}">#${r.rank}</td>`;
    html += `<td style="font-weight:600">${r.sample_id}</td>`;
    html += `<td class="score-cell">${r.score.toFixed(1)}</td>`;
    html += `<td>${h.n50 !== undefined ? h.n50.toFixed(1) : "-"}</td>`;
    html += `<td>${h.contigs !== undefined ? h.contigs.toFixed(0) : "-"}</td>`;
    html += `<td>${h.busco_pct !== undefined ? h.busco_pct.toFixed(1) + "%" : "-"}</td>`;
    html += `<td>${h.qv !== undefined ? h.qv.toFixed(1) : "-"}</td>`;
    if (hasTc) html += `<td class="${tcClass}">${tcVal !== undefined ? tcVal.toFixed(3) : "-"}</td>`;
    html += `<td>${r.best_hap}</td>`;
    html += `</tr>`;
  });

  html += `</tbody></table>`;
  $("comparison-table").innerHTML = html;
}

function selectSample(sid) { curSample = sid; renderAll(); }
function selectHap(h) { curHap = h; renderAll(); }

function metricBar(label, value, min, max, color, unit, invert) {
  if (value === undefined || value === null) return "";
  let pct = invert
    ? Math.max(0, Math.min(100, (max - value) / (max - min) * 100))
    : Math.max(0, Math.min(100, (value - min) / (max - min) * 100));
  const dv = Math.abs(value) >= 1000 ? value.toFixed(0) : (Math.abs(value) >= 100 ? value.toFixed(1) : value.toFixed(2));
  return `<div class="metric-row">
    <div class="metric-label"><span>${label}</span><span class="metric-value">${dv}${unit}</span></div>
    <div class="bar-bg"><div class="bar-fill" style="width:${pct}%;background:${color}"></div></div>
  </div>`;
}

function renderDetail(sd, ranking) {
  if (!curSample || !sd.assemblies[curSample]) { $("detail-card").innerHTML = ""; return; }

  const a = sd.assemblies[curSample];
  const ri = ranking.findIndex(r => r.sample_id === curSample);
  const rank = ri >= 0 ? ri + 1 : "?";
  const color = RANK_COLORS[Math.min(ri, RANK_COLORS.length - 1)];
  const label = ri === 0 ? "\u2605 " + RANK_LABELS[0] : RANK_LABELS[Math.min(ri, RANK_LABELS.length - 1)];
  const h = a[curHap] || {};
  const both = a.both || {};

  let html = `<div class="card selected" style="--sel-color:${color};border-color:${color}40">`;

  html += `<div class="sample-header"><div>`;
  html += `<span class="badge" style="background:${color}30;color:${color}">${label}</span>`;
  html += `<div class="sample-name">${curSample}</div></div>`;
  html += `<div class="hap-btns">`;
  ["hap1","hap2"].forEach(hp => {
    html += `<button class="hap-btn${curHap===hp?" active":""}" style="--sel-color:${color}" onclick="selectHap('${hp}')">${hp}</button>`;
  });
  html += `</div></div>`;

  html += `<div class="metrics-grid"><div>`;
  html += `<div class="metric-section-title">Contiguity</div>`;
  html += metricBar("N50", h.n50, 0, 100, color, " Mb", false);
  html += metricBar("L50 (fewer = better)", h.l50, 1, 200, color, "", true);
  html += metricBar("# Contigs (fewer = better)", h.contigs, 100, 3000, color, "", true);
  html += metricBar("Largest contig", h.largest, 0, 150, color, " Mb", false);
  html += metricBar("auN", h.aun, 0, 80, color, " Mb", false);
  html += `</div><div>`;
  html += `<div class="metric-section-title">Quality</div>`;
  html += metricBar("QV (Merqury)", h.qv, 55, 65, color, "", false);
  if (h.busco_c !== undefined && h.busco_total) {
    const bp = (h.busco_c / h.busco_total * 100);
    html += metricBar(`BUSCO complete (${bp.toFixed(1)}%)`, h.busco_c, 3100, 3450, color, "/" + h.busco_total.toFixed(0), false);
  }
  html += metricBar("BUSCO missing (fewer = better)", h.busco_m, 150, 350, color, "", true);
  html += metricBar("K-mer completeness", h.kcomp, 70, 82, color, "%", false);
  html += metricBar("HiFi depth", h.hifi_depth, 15, 45, color, "x", false);
  html += `</div></div>`;

  if (h.tc_ratio !== undefined || both.kcomp !== undefined) {
    html += `<div style="margin-top:16px"><div class="metric-section-title">Scaffolding & Combined</div>`;
    html += `<div class="metrics-grid"><div>`;
    html += metricBar("Trans/Cis ratio (lower = better)", h.tc_ratio, 0.2, 1.8, color, "", true);
    html += metricBar("Hi-C retention", h.retention_pct, 25, 40, color, "%", false);
    html += `</div><div>`;
    html += metricBar("Diploid k-mer completeness", both.kcomp, 97, 99.5, color, "%", false);
    html += metricBar("Combined QV", both.qv, 58, 65, color, "", false);
    html += `</div></div></div>`;
  }

  html += `<div class="summary-bar">`;
  html += `<span class="hl" style="color:${color}">Total length:</span> ${h.total !== undefined ? (h.total/1000).toFixed(2) + " Gb" : "-"} &nbsp;|&nbsp;`;
  html += `<span class="hl" style="color:${color}">GC:</span> ${h.gc !== undefined ? h.gc.toFixed(2) + "%" : "-"} &nbsp;|&nbsp;`;
  html += `<span class="hl" style="color:${color}">N's/100kb:</span> ${h.ns_100k !== undefined ? h.ns_100k.toFixed(2) : "-"} &nbsp;|&nbsp;`;
  html += `<span class="hl" style="color:${color}">BUSCO dup:</span> ${h.busco_d !== undefined ? h.busco_d.toFixed(0) : "-"}`;
  html += `</div>`;

  html += `</div>`;
  $("detail-card").innerHTML = html;
}

init();
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Generate assembly QC HTML report")
    parser.add_argument("--input",  required=True, help="Input metrics CSV")
    parser.add_argument("--output", required=True, help="Output HTML file")
    parser.add_argument("--report-stage", default="final",
                        help="Default stage to display (default: final)")
    args = parser.parse_args()

    data = parse_metrics(args.input)
    if not data:
        print("ERROR: No metrics parsed from input CSV", file=sys.stderr)
        sys.exit(1)

    stages = get_available_stages(data)
    if not stages:
        print("ERROR: No valid stages found", file=sys.stderr)
        sys.exit(1)

    print(f"Parsed {len(data)} samples across {len(stages)} stages: {', '.join(stages)}",
          file=sys.stderr)

    report = build_report_data(data, stages)

    report_json = json.dumps(report, indent=None, separators=(",", ":"))
    html = HTML_TEMPLATE.replace("__REPORT_JSON__", report_json)
    html = html.replace("__DEFAULT_STAGE__", args.report_stage)

    with open(args.output, "w") as fh:
        fh.write(html)

    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()