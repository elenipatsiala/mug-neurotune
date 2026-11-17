// µG-NeuroTune scoring core (browser + node). Must match mug_neurotune/scoring.py exactly.
(function (root) {
  function parseTable(text) {
    const lines = text.replace(/\r/g, "").split("\n").filter(l => l.trim().length && !l.startsWith("#"));
    const sep = (lines[0].match(/\t/g) || []).length > (lines[0].match(/,/g) || []).length ? "\t" : ",";
    const head = lines[0].split(sep).map(s => s.trim().replace(/^"|"$/g, ""));
    const rows = lines.slice(1).map(l => l.split(sep).map(s => s.trim().replace(/^"|"$/g, "")));
    return { head, rows };
  }

  // counts: genes x samples. Duplicate symbols summed. Non-numeric annotation columns dropped.
  function parseCounts(text) {
    const { head, rows } = parseTable(text);
    const skip = new Set(["Chr", "Start", "End", "Strand", "Length", "gene_name", "Symbol", "symbol"]);
    let symCol = head.findIndex(h => ["gene_name", "Symbol", "symbol"].includes(h));
    const sampleIdx = [];
    for (let j = 1; j < head.length; j++) if (!skip.has(head[j])) sampleIdx.push(j);
    const samples = sampleIdx.map(j => head[j]);
    const map = new Map();
    for (const r of rows) {
      const g = String(symCol > 0 ? r[symCol] : r[0]);
      if (!g) continue;
      const v = sampleIdx.map(j => { const x = parseFloat(r[j]); return isFinite(x) ? x : 0; });
      if (map.has(g)) { const a = map.get(g); for (let k = 0; k < v.length; k++) a[k] += v[k]; }
      else map.set(g, v);
    }
    const genes = [...map.keys()].sort();
    return { genes, samples, X: genes.map(g => map.get(g)) };
  }

  function parseMeta(text) {
    const { head, rows } = parseTable(text);
    const out = {};
    for (const r of rows) { const o = {}; head.forEach((h, j) => o[h] = r[j]); out[r[0]] = o; }
    return { columns: head, bySample: out };
  }

  function logCPM(c) {
    const n = c.samples.length, lib = new Array(n).fill(0);
    for (const row of c.X) for (let j = 0; j < n; j++) lib[j] += row[j];
    return c.X.map(row => row.map((v, j) => Math.log2(v / lib[j] * 1e6 + 1)));
  }

  // average-tie rank / n, per sample (column)
  function rankColumns(L) {
    const G = L.length, n = L[0].length, R = Array.from({ length: G }, () => new Array(n));
    for (let j = 0; j < n; j++) {
      const idx = [...Array(G).keys()].sort((a, b) => L[a][j] - L[b][j]);
      let i = 0;
      while (i < G) {
        let k = i;
        while (k + 1 < G && L[idx[k + 1]][j] === L[idx[i]][j]) k++;
        const avg = (i + k) / 2 + 1;
        for (let t = i; t <= k; t++) R[idx[t]][j] = avg / G;
        i = k + 1;
      }
    }
    return R;
  }

  const lowerBound = (a, x) => { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < x) lo = m + 1; else hi = m; } return lo; };
  const upperBound = (a, x) => { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] <= x) lo = m + 1; else hi = m; } return lo; };

  // Background-referenced rank score (matches scoring.py): each program gene is ranked
  // against genes that belong to NO program, so one program rising cannot drag others down.
  function scoreAll(c, genesets, clock) {
    const L = logCPM(c), gi = new Map(c.genes.map((g, i) => [g, i]));
    const programGenes = new Set(Object.values(genesets).flatMap(d => d.genes));
    let ckUp = [], ckDown = [];
    if (clock) {
      const res = g => gi.has(g) ? g : (clock.aliases[g] && gi.has(clock.aliases[g]) ? clock.aliases[g] : null);
      ckUp = clock.up.map(res).filter(Boolean); ckDown = clock.down.map(res).filter(Boolean);
    }
    let bgIdx = c.genes.map((g, i) => programGenes.has(g) ? -1 : i).filter(i => i >= 0);
    if (bgIdx.length < 100) bgIdx = c.genes.map((_, i) => i);
    const n = c.samples.length;
    const bg = [...Array(n).keys()].map(j => bgIdx.map(i => L[i][j]).sort((a, b) => a - b));
    const scores = {}, coverage = {};
    for (const [name, def] of Object.entries(genesets)) {
      const present = def.genes.filter(g => gi.has(g));
      coverage[name] = { n: def.genes.length, present: present.length, missing: def.genes.filter(g => !gi.has(g)) };
      scores[name] = c.samples.map((_, j) => {
        if (!present.length) return NaN;
        const ref = bg[j]; let s = 0;
        for (const g of present) { const v = L[gi.get(g)][j]; s += (lowerBound(ref, v) + upperBound(ref, v)) / 2 / ref.length; }
        return s / present.length - 0.5;
      });
    }
    if (scores.da_identity && scores.off_target)
      scores.fidelity = scores.da_identity.map((v, j) => v - scores.off_target[j]);
    if (clock) {
      coverage.age_clock = { n: clock.up.length + clock.down.length, present: ckUp.length + ckDown.length, missing: [] };
      if (ckUp.length >= 0.8 * clock.up.length && ckDown.length >= 0.8 * clock.down.length) {
        const ex = new Set([...programGenes, ...ckUp, ...ckDown]);
        const bIdx = c.genes.map((g, i) => ex.has(g) ? -1 : i).filter(i => i >= 0);
        const bg2 = [...Array(n).keys()].map(j => bIdx.map(i => L[i][j]).sort((a, b) => a - b));
        const mr = (genes, j) => genes.reduce((s2, g) => { const v = L[gi.get(g)][j], r = bg2[j]; return s2 + (lowerBound(r, v) + upperBound(r, v)) / 2 / r.length; }, 0) / genes.length;
        scores.age_clock = c.samples.map((_, j) => mr(ckUp, j) - mr(ckDown, j));
      }
    }
    // validated progenitor-loss read-out: minus the program-style rank score of SOX2 + NES
    if (gi.has("SOX2") && gi.has("NES")) {
      const ref = bg;
      scores.progenitor_loss = c.samples.map((_, j) => {
        let s2 = 0; for (const g of ["SOX2", "NES"]) { const v = L[gi.get(g)][j]; s2 += (lowerBound(ref[j], v) + upperBound(ref[j], v)) / 2 / ref[j].length; }
        return -(s2 / 2 - 0.5);
      });
    }
    return { scores, coverage };
  }

  // Within-stratum (line|microglia) effects, test - ref, with a within-stratum permutation
  // test and Benjamini-Hochberg q-values across programs (same logic as compare.py).
  function stratifiedEffects(samples, scores, meta, ref, test, nPerm = 2000, seed = 1) {
    const strata = {};
    samples.forEach((s, j) => {
      const m = meta.bySample[s]; if (!m || (m.condition !== ref && m.condition !== test)) return;
      const key = [m.line, m.microglia].join("|");
      (strata[key] = strata[key] || []).push({ j, c: m.condition });
    });
    const keys = Object.keys(strata), groups = keys.map(k => strata[k]);
    const lineOf = keys.map(k => k.split("|")[0]);
    let st = seed >>> 0; const rnd = () => ((st = (st * 1664525 + 1013904223) >>> 0) / 4294967296);
    const shuffle = a => { a = a.slice(); for (let i = a.length - 1; i > 0; i--) { const k = Math.floor(rnd() * (i + 1)); [a[i], a[k]] = [a[k], a[i]]; } return a; };
    function eff(vals, labels) {
      let dsum = 0, w = 0, gsum = 0, wg = 0;
      groups.forEach((members, gix) => {
        const lab = labels[gix];
        const a = members.filter((_, t) => lab[t] === ref).map(x => vals[x.j]);
        const b = members.filter((_, t) => lab[t] === test).map(x => vals[x.j]);
        if (!a.length || !b.length) return;
        const ma = mean(a), mb = mean(b), n = a.length + b.length;
        dsum += (mb - ma) * n; w += n;
        if (a.length > 1 && b.length > 1) {
          const sp = Math.sqrt(((a.length - 1) * vr(a) + (b.length - 1) * vr(b)) / (n - 2));
          if (sp > 0) { const J = 1 - 3 / (4 * n - 9); gsum += J * (mb - ma) / sp * n; wg += n; }
        }
      });
      return { diff: w ? dsum / w : NaN, g: wg ? gsum / wg : NaN };
    }
    const obsLabels = groups.map(m => m.map(x => x.c));
    const permLabels = [...Array(nPerm)].map(() => obsLabels.map(shuffle));
    const out = {};
    for (const [prog, vals] of Object.entries(scores)) {
      const o = eff(vals, obsLabels); let ext = 0, tot = 0;
      for (const pl of permLabels) { const e = eff(vals, pl).g; if (isFinite(e)) { tot++; if (Math.abs(e) >= Math.abs(o.g) - 1e-12) ext++; } }
      // direction must agree in every cell line
      const perLine = {};
      groups.forEach((members, gix) => {
        const a = members.filter(x => x.c === ref).map(x => vals[x.j]);
        const b = members.filter(x => x.c === test).map(x => vals[x.j]);
        if (a.length && b.length) (perLine[lineOf[gix]] = perLine[lineOf[gix]] || []).push(mean(b) - mean(a));
      });
      const consistent = Object.values(perLine).every(d => Math.sign(mean(d)) === Math.sign(o.diff));
      out[prog] = { ...o, p: (ext + 1) / (tot + 1), consistent };
    }
    const progs = Object.keys(out).sort((a, b) => out[b].p - out[a].p); let prev = 1;
    progs.forEach((k, r) => { const rank = progs.length - r; prev = Math.min(prev, out[k].p * progs.length / rank); out[k].q = prev; });
    return { effects: out, nStrata: groups.length };
  }
  const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
  const vr = a => { const m = mean(a); return a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1); };

  const api = { parseCounts, parseMeta, logCPM, rankColumns, scoreAll, stratifiedEffects };
  if (typeof module !== "undefined") module.exports = api; else root.MUGNT = api;
})(this);
