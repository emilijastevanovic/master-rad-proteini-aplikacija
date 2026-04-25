let SEQ_CACHE = null;
let SEQ_RENDERED = false;

function showViewerMessage(msg, isError = false) {
  const el = document.getElementById("viewer-message");
  if (!el) return;
  el.textContent = msg;
  el.className = isError ? "error" : "";
  el.style.visibility = msg ? "visible" : "hidden";
}

function hideViewerMessage() {
  showViewerMessage("");
}

let HEATMAP_URLS = null;
let SEGMENTS_URL = null;

let TOOLTIP_LOCKED = false;

let LAST_ANALYSIS_KEY = null;
function analysisKey(protein, chain, type, maxdist, mindist) {
  return [protein, chain, type, maxdist, mindist].join("|");
}

// =========================
//   SHOW TOOLTIP FUNCTION
// =========================


function showSphereTooltip(pixelX, pixelY, html) {
  const tooltip = document.getElementById("nodeTooltip");
  const wrapper = document.getElementById("viewer3d-wrapper");

  if (!tooltip) return;

  tooltip.innerHTML = html;

  // ako nemamo validne koordinate klika, prikaži tooltip u okviru viewera
  if (!Number.isFinite(pixelX) || !Number.isFinite(pixelY) || (pixelX === 0 && pixelY === 0)) {
    if (wrapper) {
      const r = wrapper.getBoundingClientRect();
      tooltip.style.left = (r.left + 20) + "px";
      tooltip.style.top  = (r.top + 20) + "px";
    } else {
      tooltip.style.left = "20px";
      tooltip.style.top  = "20px";
    }
  } else {
    tooltip.style.left = (pixelX + 12) + "px";
    tooltip.style.top  = (pixelY + 12) + "px";
  }

  tooltip.style.visibility = "visible";
}


// =========================
//   3D GRAPH
// =========================
async function loadGraph3D() {
  console.log("loadGraph3D CALLED");
  const tooltip = document.getElementById("nodeTooltip");
  if (tooltip) tooltip.style.visibility = "hidden";

  const loader = document.getElementById("loader3d");
  if (loader) loader.style.visibility = "visible";

  const el = document.getElementById("viewer3d");
  if (el) el.innerHTML = "";

  const protein   = document.getElementById("protein")?.value.trim() || "";
  const aminoname = document.getElementById("aminoname")?.value.trim() || "";
  const type      = document.getElementById("type")?.value || "prazno";
  const maxdist   = document.getElementById("maxdist")?.value.trim() || "";
  const mindist   = document.getElementById("mindist")?.value.trim() || "";
  const chain     = document.getElementById("chain")?.value.trim() || "";

  if (!protein) {
    showViewerMessage("Unesi naziv proteina.", true);
    if (loader) loader.style.visibility = "hidden";
    return;
  }
  hideViewerMessage();

  const key = analysisKey(protein, chain, type, maxdist, mindist);
  if (key !== LAST_ANALYSIS_KEY) {
    SEQ_CACHE = null;
    SEQ_RENDERED = false;
    HEATMAP_URLS = null;
    SEGMENTS_URL = null;
    LAST_ANALYSIS_KEY = key;

    // očisti vidljiv sadržaj tabova
    document.querySelectorAll(".aa-missing").forEach(x => x.remove());
    const distTable = document.getElementById("distSummaryTable");
    if (distTable) distTable.style.display = "none";
    const distMeta = document.getElementById("distSummaryMeta");
    if (distMeta) distMeta.textContent = "Klikni \"Prikaži\" da učitaš podatke.";
    const seqEl = document.getElementById("aaSequenceTokens");
    const ssEl  = document.getElementById("ssTrackTokens");
    if (seqEl) seqEl.textContent = "Učitavanje...";
    if (ssEl)  ssEl.textContent  = "";

    // reset Outliers tab
    const outlierSSMeta  = document.getElementById("outlierSSMeta");
    const outlierTopMeta = document.getElementById("outlierTopMeta");
    const outlierSSTable  = document.getElementById("outlierSSTable");
    const outlierTopTable = document.getElementById("outlierTopTable");
    if (outlierSSMeta)  outlierSSMeta.textContent  = "Klikni \"Prikaži\" da učitaš podatke.";
    if (outlierTopMeta) outlierTopMeta.textContent = "";
    if (outlierSSTable)  { outlierSSTable.style.display  = "none"; outlierSSTable.querySelector("tbody").innerHTML  = ""; }
    if (outlierTopTable) { outlierTopTable.style.display = "none"; outlierTopTable.querySelector("tbody").innerHTML = ""; }

    document.querySelectorAll("img.heatmap-img, img.seg-img, .img-loading").forEach(x => x.remove());
    const heatCard = document.querySelector("#tabDistance .card");
    const segCard  = document.querySelector("#tabSegments .card");
    if (heatCard && !heatCard.querySelector(".muted")) {
      const p = document.createElement("div");
      p.className = "muted";
      p.textContent = "Klikni \"Prikaži\" pa otvori ovaj tab.";
      heatCard.appendChild(p);
    }
    if (segCard && !segCard.querySelector(".muted")) {
      const p = document.createElement("div");
      p.className = "muted";
      p.textContent = "Klikni \"Prikaži\" pa otvori ovaj tab.";
      segCard.appendChild(p);
    }
  }

  // STATISTIKE
  loadAAComposition(protein, chain, type, maxdist, mindist);
  loadSSPairs(protein, chain, type, maxdist, mindist);
  loadDistanceSummary(protein, chain, type, maxdist, mindist);
  loadOutlierSS(protein, chain, type, maxdist, mindist);

  prefetchSeq(protein, chain, type, maxdist, mindist);

  setHeatmapUrls(protein, chain, type, maxdist, mindist);
  setSegmentsUrl(protein, chain);

  // ako je tab već otvoren, odmah re-renderuj
  if (document.getElementById("tabDistance")?.style.display !== "none") {
    renderHeatmapsFromUrls();
  }
  if (document.getElementById("tabSegments")?.style.display !== "none") {
    renderSegmentsImg();
  }
  if (document.getElementById("tabSeq")?.style.display !== "none") {
    renderSeqIfReady();
  }

  try {
    // ✅ važna stvar: encodeURIComponent za svaki parametar
    const url =
      `/api/graph3d/?protein=${encodeURIComponent(protein)}` +
      `&aminoname=${encodeURIComponent(aminoname)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&max_distance=${encodeURIComponent(maxdist)}` +
      `&min_distance=${encodeURIComponent(mindist)}`;

    console.log("Fetching graph:", url);

    const res = await fetch(url);
    if (!res.ok) {
      throw new Error(`HTTP ${res.status} ${res.statusText}`);
    }

    const json = await res.json();
    const nodes = json.nodes || [];
    const edges = json.edges || [];

    // Ako nema viewer elementa, nema smisla dalje
    if (!el) {
      console.warn("viewer3d element not found. Cannot render 3D.");
      return;
    }
    if (!nodes.length) {
      showViewerMessage("Nema rezultata za izabrane filtere.");
      return;
    }

    const viewer = $3Dmol.createViewer(el, { backgroundColor: "white" });

    // zaštita ako nema koordinata
    const xs = nodes.map(n => n.x).filter(v => typeof v === "number");
    const ys = nodes.map(n => n.y).filter(v => typeof v === "number");
    const zs = nodes.map(n => n.z).filter(v => typeof v === "number");

    if (!xs.length || !ys.length || !zs.length) {
      console.warn("No coordinates in nodes - cannot render 3D.");
      viewer.render();
      return;
    }

    const min = a => Math.min(...a);
    const max = a => Math.max(...a);

    const cx = (min(xs) + max(xs)) / 2;
    const cy = (min(ys) + max(ys)) / 2;
    const cz = (min(zs) + max(zs)) / 2;

    const span = Math.max(max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)) || 1;
    const SCALE = 60 / span;

    const nodeById = {};
    nodes.forEach(n => nodeById[n.id] = n);

    const aaColorMap3d = {
      // hidrofobne — žuta
      ALA: 0xe6b800, VAL: 0xe6b800, ILE: 0xe6b800, LEU: 0xe6b800,
      MET: 0xe6b800, PHE: 0xe6b800, TRP: 0xe6b800, PRO: 0xe6b800,
      // polarne — zelena
      SER: 0x4caf50, THR: 0x4caf50, CYS: 0x4caf50, TYR: 0x4caf50,
      ASN: 0x4caf50, GLN: 0x4caf50,
      // pozitivne — plava
      LYS: 0x2979ff, ARG: 0x2979ff, HIS: 0x2979ff,
      // negativne — crvena
      ASP: 0xe53935, GLU: 0xe53935,
      // GLY — siva
      GLY: 0x9e9e9e,
    };
    const defaultColor = 0x9e9e9e;

    // spheres
    nodes.forEach(n => {
      n._x = (n.x - cx) * SCALE;
      n._y = (n.y - cy) * SCALE;
      n._z = (n.z - cz) * SCALE;

      n._aa = n.aa;
      n._index = n.index;
      n._protein = n.protein;
      n._chain = n.chain;
      n._ss = n.ss;

      viewer.addSphere({
        center: { x: n._x, y: n._y, z: n._z },
        radius: 0.7,
        color: aaColorMap3d[n.aa] || defaultColor,
        clickable: true,
        callback: function(atom, viewer, event) {
          console.log("SPHERE CALLBACK FIRED", n._aa, n._index, event);
          if (event) event.stopPropagation();
          
          TOOLTIP_LOCKED = true;
          const x = (event && typeof event.clientX === "number") ? event.clientX : NaN;
          const y = (event && typeof event.clientY === "number") ? event.clientY : NaN;
        
          const html = `
            <b>${n._aa}${n._index}</b><br>
            Protein: ${n._protein}<br>
            Lanac: ${n._chain}<br>
            SS: ${n._ss || "-"}
          `;
          showSphereTooltip(x, y, html);
        }
      });
    });

    // edges
    edges.forEach(e => {
      const a = nodeById[e.source];
      const b = nodeById[e.target];
      if (!a || !b) return;

      const edgeColor =
        e.type === "caca"    ? 0x1f77b4 :
        e.type === "minbezh" ? 0x2ca02c :
        e.type === "maxbezh" ? 0xd62728 :
                               0x444444;

      viewer.addCylinder({
        start:    { x: a._x, y: a._y, z: a._z },
        end:      { x: b._x, y: b._y, z: b._z },
        radius:   0.08,
        color:    edgeColor,
        clickable: true,
        callback: function(shape, viewer, event) {
          if (event) event.stopPropagation();
          TOOLTIP_LOCKED = true;
          const x = (event && typeof event.clientX === "number") ? event.clientX : NaN;
          const y = (event && typeof event.clientY === "number") ? event.clientY : NaN;

          const html = `
            <b>${a._aa}${a._index} — ${b._aa}${b._index}</b><br>
            Tip: ${e.type || "–"}<br>
            Rastojanje: ${e.value != null ? Number(e.value).toFixed(2) + " Å" : "–"}<br>
            Lanac: ${a._chain}
          `;
          showSphereTooltip(x, y, html);
        }
      });
    });

    viewer.zoomTo();
    viewer.render();

  } catch (err) {
    console.error("loadGraph3D failed", err);
    showViewerMessage("Greška pri učitavanju grafa.", true);
  } finally {
    if (loader) loader.style.visibility = "hidden";
  }
}


async function loadAAComposition(protein, chain, type, maxdist, mindist) {
  const meta = document.getElementById("aaStatsMeta");
  const tbody = document.querySelector("#aaStatsTable tbody");

  if (!meta || !tbody) return;

  // reset UI
  meta.textContent = "Učitavanje...";
  tbody.innerHTML = "";

  try {
    const url = `/api/stats/?protein=${encodeURIComponent(protein)}&type=${encodeURIComponent(type)}&include=aa_composition,sequence`;
    const res = await fetch(url);

    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();

    // ✅ očekujemo: data.aa_composition = { protein, total, rows }
    const block = data.aa_composition;
    const rows = block.aa_composition;
    const total = block.total || 0;

    if (!rows.length) {
      meta.textContent = "Nema podataka za izabrani protein.";
      return;
    }

  
    for (const r of rows) {
      const aa = r.name ?? r.aa ?? "";
      const count = Number(r.cnt ?? r.count ?? 0);

      const pct =  (total > 0 ? (count / total) * 100 : 0);

      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${aa}</td>
        <td class="num">${count}</td>
        <td class="num">${pct.toFixed(1)}%</td>
      `;
      tbody.appendChild(tr);
    }

    const missing = data.aa_composition?.missing || [];
    if (missing.length) {
      const existing = tbody.closest("table").nextElementSibling;
      if (existing?.classList.contains("aa-missing")) existing.remove();

      const missingEl = document.createElement("div");
      missingEl.className = "muted aa-missing";
      missingEl.style.marginTop = "8px";
      missingEl.style.fontSize = "0.85rem";
      missingEl.textContent = `Nije prisutna: ${missing.join(", ")}`;
      tbody.closest("table").insertAdjacentElement("afterend", missingEl);
    }

    meta.textContent = `Ukupno aminokiselina: ${total} · Prisutnih tipova: ${rows.length} / 20`;
  } catch (err) {
    console.error("loadAAComposition failed ❌", err);
    meta.textContent = "Greška pri učitavanju (proveri konzolu).";
  }
}

async function prefetchSeq(protein, chain, type, maxdist, mindist) {


  if (SEQ_CACHE ) return;

  SEQ_CACHE = null;
  SEQ_RENDERED = false;

  try {
    const url = `/api/stats/?protein=${encodeURIComponent(protein)}&type=${encodeURIComponent(type)}&include=aa_composition,sequence`;
    const res = await fetch(url);

    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();

    // očekujemo data.sequence = { protein, aa_seq, aa_ss_seq }
    const block = data.sequence;
    if (!block) return;

    // **OVDE SPajamo AA + SS u jednu strukturu po poziciji**
    // ako backend vraća liste rows1/rows2, onda:
    //   aa_seq = [{name:"PRO"}, ...]
    //   aa_ss_seq = [{ss:"H"}, ...]
    // a ako vrati string, onda mora drugačije (spomenuću dole)
    const rows1 = block.aa_seq || [];
    const rows2 = block.aa_ss_seq || [];

    const n = Math.min(rows1.length, rows2.length);

    SEQ_CACHE = Array.from({ length: n }, (_, i) => ({
      aa: rows1[i]?.name ?? "",
      ss: rows2[i]?.ss ?? "."
    }));

  } catch (err) {
    console.error("prefetchSeq failed", err);
    SEQ_CACHE = null;
  }
}

function renderSeqIfReady() {
  const seqEl = document.getElementById("aaSequenceTokens");
  const ssEl = document.getElementById("ssTrackTokens");

  if (!seqEl) return;

  // ako još ništa nije prefetched
  if (!SEQ_CACHE) {
    seqEl.textContent = "Klikni “Prikaži” da učitaš podatke.";
    if (ssEl) ssEl.textContent = "";
    return;
  }

  // ako je već renderovano, samo izađi
  if (SEQ_RENDERED) return;

  // render AA tokens
  seqEl.innerHTML = "";
  for (const item of SEQ_CACHE) {
    const aa = item.aa || "---";
    const ss = item.ss || ".";

    const span = document.createElement("span");
    span.className = `seqToken ${aaToClass(aa)}`;
    span.textContent = `[${aa}]`;
    span.title = `AA=${aa}, SS=${ss}`;
    seqEl.appendChild(span);
  }

  // render SS track (opciono, kao [H] [E] ...)
  if (ssEl) {
    ssEl.innerHTML = "";
    for (const item of SEQ_CACHE) {
      const ss = item.ss || ".";

      const span = document.createElement("span");
      span.className = `seqToken ${ssToClass(ss)}`;
      span.textContent = `[${ss}]`;
      span.title = `SS=${ss}`;
      ssEl.appendChild(span);
    }
  }

  SEQ_RENDERED = true;
}

function setHeatmapUrls(protein, chain, type, maxdist, mindist) {
  const p = encodeURIComponent(protein);
  const c = encodeURIComponent(chain);
  const t = encodeURIComponent(type);
  const maxd = encodeURIComponent(maxdist);
  const mind = encodeURIComponent(mindist);

  const stamp = Date.now(); // cache-bust

  HEATMAP_URLS = {
    distance: `/api/heatmap/distance/?protein=${p}&chain=${c}&type=${t}&max_distance=${maxd}&min_distance=${mind}&_=${stamp}`,
    meta: { protein, chain, type, maxdist, mindist }
  };
}

function setSegmentsUrl(protein, chain) {
  const stamp = Date.now();
  SEGMENTS_URL = `/api/heatmap/segments/?protein=${encodeURIComponent(protein)}&chain=${encodeURIComponent(chain)}&_=${stamp}`;
}

function renderSegmentsImg() {
  const tab = document.getElementById("tabSegments");
  if (!tab) return;
  const card = tab.querySelector(".card");
  if (!card) return;

  if (!SEGMENTS_URL) {
    const muted = card.querySelector(".muted");
    if (muted) muted.textContent = "Klikni \"Prikaži\" da učitaš histogram.";
    return;
  }

  const muted = card.querySelector(".muted");
  if (muted) muted.remove();
  card.querySelectorAll("img.seg-img").forEach(x => x.remove());

  const loadingSeg = document.createElement("div");
  loadingSeg.className = "img-loading";
  loadingSeg.textContent = "Učitavanje...";
  card.appendChild(loadingSeg);

  const img = document.createElement("img");
  img.className = "seg-img zoomable";
  img.src = SEGMENTS_URL;
  img.alt = "Segment length histogram";
  img.style.width = "100%";
  img.style.borderRadius = "10px";
  img.style.marginTop = "10px";
  img.onload  = () => loadingSeg.remove();
  img.onerror = () => {
    loadingSeg.remove();
    img.remove();
    const err = document.createElement("div");
    err.className = "muted";
    err.textContent = "Nema podataka za segmente (klikni Prikaži).";
    card.appendChild(err);
  };
  card.appendChild(img);
}

function renderHeatmapsFromUrls() {
  const tab = document.getElementById("tabDistance");
  if (!tab) return;

  const card = tab.querySelector(".card");
  if (!card) return;

  if (!HEATMAP_URLS) {
    const muted = card.querySelector(".muted");
    if (muted) muted.textContent = "Klikni “Prikaži” da učitaš heatmap.";
    return;
  }

  // ukloni placeholder
  const muted = card.querySelector(".muted");
  if (muted) muted.remove();

  // obriši prethodne slike
  card.querySelectorAll("img.heatmap-img").forEach(x => x.remove());
  card.querySelectorAll(".heatmap-subtitle").forEach(x => x.remove());

  const { distance } = HEATMAP_URLS;

  const loading1 = document.createElement("div");
  loading1.className = "img-loading";
  loading1.textContent = "Učitavanje...";
  card.appendChild(loading1);

  const img1 = document.createElement("img");
  img1.className = "heatmap-img zoomable";
  img1.src = distance;
  img1.alt = "Distance heatmap";
  img1.style.width = "100%";
  img1.style.borderRadius = "10px";
  img1.style.marginTop = "10px";
  img1.onload  = () => loading1.remove();
  img1.onerror = () => {
    loading1.remove();
    img1.remove();
    const err = document.createElement("div");
    err.className = "muted";
    err.textContent = "Nema podataka za heatmap (izaberi tip rastojanja i klikni Prikaži).";
    card.appendChild(err);
  };
  card.appendChild(img1);
}


async function loadDistanceSummary(protein, chain, type, maxdist, mindist) {
  const meta  = document.getElementById("distSummaryMeta");
  const table = document.getElementById("distSummaryTable");
  const tbody = table?.querySelector("tbody");

  if (!meta || !table || !tbody) return;

  meta.textContent = "Učitavanje...";
  table.style.display = "none";
  tbody.innerHTML = "";

  try {
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&max_distance=${encodeURIComponent(maxdist)}` +
      `&min_distance=${encodeURIComponent(mindist)}` +
      `&include=distance_summary`;

    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();
    const s = data.distance_summary;

    if (!s || s.empty) {
      meta.textContent = "Nema podataka za izabrani tip.";
      return;
    }

    const fmt = v => (v == null ? "–" : Number(v).toFixed(3));
    const rows = [
      ["Tip",       s.type],
      ["Count",     s.count],
      ["Mean",      fmt(s.mean) + " Å"],
      ["Median",    fmt(s.median) + " Å"],
      ["p95",       fmt(s.p95) + " Å"],
      ["Min",       fmt(s.min) + " Å"],
      ["Max",       fmt(s.max) + " Å"],
      ["Pct < 8Å",  s.pct_lt_8 != null ? s.pct_lt_8.toFixed(1) + " %" : "–"],
    ];

    for (const [label, value] of rows) {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td>${label}</td><td class="num">${value}</td>`;
      tbody.appendChild(tr);
    }

    meta.textContent = "";
    table.style.display = "";

  } catch (err) {
    console.error("loadDistanceSummary failed", err);
    meta.textContent = "Greška pri učitavanju (proveri konzolu).";
  }
}

async function loadSSPairs(protein, chain, type, maxdist, mindist) {
  const meta  = document.getElementById("ssPairsMeta");
  const table = document.getElementById("ssPairsTable");
  const tbody = table?.querySelector("tbody");

  if (!meta || !table || !tbody) return;

  meta.textContent = "Učitavanje...";
  table.style.display = "none";
  tbody.innerHTML = "";

  try {
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&include=ss_pairs`;

    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();
    const rows = data.ss_pairs?.ss_pairs || [];

    if (!rows.length) {
      meta.textContent = "Nema podataka za izabrani protein/lanac/tip.";
      return;
    }

    const fmt = v => (v == null ? "–" : Number(v).toFixed(2));

    for (const r of rows) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${r.ss_pair}</td>
        <td class="num">${r.count}</td>
        <td class="num">${fmt(r.mean)}</td>
        <td class="num">${fmt(r.median)}</td>
        <td class="num">${fmt(r.min)}</td>
        <td class="num">${fmt(r.max)}</td>
        <td class="num">${fmt(r.std)}</td>
      `;
      tbody.appendChild(tr);
    }

    meta.textContent = `Tip: ${data.ss_pairs?.type ?? "–"} · Parova: ${rows.length}`;
    table.style.display = "";

  } catch (err) {
    console.error("loadSSPairs failed", err);
    meta.textContent = "Greška pri učitavanju (proveri konzolu).";
  }
}

async function loadOutlierSS(protein, chain, type, maxdist, mindist) {
  const ssMeta   = document.getElementById("outlierSSMeta");
  const ssTable  = document.getElementById("outlierSSTable");
  const ssTbody  = ssTable?.querySelector("tbody");
  const topMeta  = document.getElementById("outlierTopMeta");
  const topTable = document.getElementById("outlierTopTable");
  const topTbody = topTable?.querySelector("tbody");

  if (!ssMeta || !ssTable || !ssTbody || !topMeta || !topTable || !topTbody) return;

  ssMeta.textContent  = "Učitavanje...";
  topMeta.textContent = "";
  ssTable.style.display  = "none";
  topTable.style.display = "none";
  ssTbody.innerHTML  = "";
  topTbody.innerHTML = "";

  try {
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&include=outlier_ss`;

    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();
    const block = data.outlier_ss;

    const ssStats     = block?.ss_stats     || [];
    const topResidues = block?.top_residues || [];

    const fmt = v => (v == null ? "–" : Number(v).toFixed(3));

    if (!ssStats.length) {
      ssMeta.textContent = "Nema podataka za izabrani protein/lanac/tip.";
    } else {
      for (const r of ssStats) {
        const tr = document.createElement("tr");
        tr.innerHTML = `
          <td>${r.ss ?? "–"}</td>
          <td class="num">${r.n_rezidua ?? "–"}</td>
          <td class="num">${fmt(r.mean)}</td>
          <td class="num">${fmt(r.median)}</td>
          <td class="num">${fmt(r.p90)}</td>
          <td class="num">${fmt(r.max)}</td>
        `;
        ssTbody.appendChild(tr);
      }
      ssMeta.textContent = `Tip: ${block?.type ?? "–"} · SS tipova: ${ssStats.length}`;
      ssTable.style.display = "";
    }

    if (topResidues.length) {
      for (const r of topResidues) {
        const tr = document.createElement("tr");
        tr.innerHTML = `
          <td>${r.aa ?? "–"}</td>
          <td class="num">${r.index ?? "–"}</td>
          <td>${r.chain ?? "–"}</td>
          <td>${r.ss ?? "–"}</td>
          <td class="num">${fmt(r.avg_dist)}</td>
        `;
        topTbody.appendChild(tr);
      }
      topMeta.textContent = `Prikazano top ${topResidues.length} najizolovanijih.`;
      topTable.style.display = "";
    }

  } catch (err) {
    console.error("loadOutlierSS failed", err);
    ssMeta.textContent = "Greška pri učitavanju (proveri konzolu).";
  }
}

function ssToClass(ss) {
  if (ss === "H") return "ss-H";
  if (ss === "E") return "ss-E";
  if (ss === "T") return "ss-T";
  if (ss === "S") return "ss-S";
  return "ss-dot";
}

const AA_HYDROPHOBIC = new Set(["ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO"]);
const AA_POLAR       = new Set(["SER","THR","CYS","TYR","ASN","GLN"]);
const AA_POSITIVE    = new Set(["LYS","ARG","HIS"]);
const AA_NEGATIVE    = new Set(["ASP","GLU"]);

function aaToClass(aa) {
  if (AA_HYDROPHOBIC.has(aa)) return "aa-hydrophobic";
  if (AA_POLAR.has(aa))       return "aa-polar";
  if (AA_POSITIVE.has(aa))    return "aa-positive";
  if (AA_NEGATIVE.has(aa))    return "aa-negative";
  return "aa-special";
}

function setupTabs() {
  const buttons = document.querySelectorAll(".tab-btn");
  const tabs = document.querySelectorAll(".tab-content");

  buttons.forEach(btn => {
    btn.addEventListener("click", () => {
      const tabId = btn.getAttribute("data-tab");
      if (!tabId) return;

      // 1) active klasa na dugmad
      buttons.forEach(b => b.classList.remove("active"));
      btn.classList.add("active");

      // 2) sakrij sve tabove
      tabs.forEach(t => (t.style.display = "none"));

      // 3) prikaži izabrani
      const activeTab = document.getElementById(tabId);
      if (activeTab) activeTab.style.display = "block";

      // 4) ako je Seq tab, renderuj kad se prvi put otvori
      if (tabId === "tabSeq") {
        renderSeqIfReady();
      }
      if (tabId === "tabDistance") {
        renderHeatmapsFromUrls();
      }
      if (tabId === "tabSegments") {
        renderSegmentsImg();
      }
    });
  });
}

function hideNodeTooltip() {
  const tooltip = document.getElementById("nodeTooltip");
  if (tooltip) tooltip.style.visibility = "hidden";
}

// Pozovi setupTabs kad se DOM učita
document.addEventListener("DOMContentLoaded", () => {
  // 1) inicijalizuj tabove (ovo ti je falilo)
  setupTabs();

  // 2) tooltip hide logika
  const tooltip = document.getElementById("nodeTooltip");
  const wrapper = document.getElementById("viewer3d-wrapper");
  const viewer  = document.getElementById("viewer3d");

  if (tooltip) {
    tooltip.style.pointerEvents = "none";
  }

  // klik van wrapper-a
  document.addEventListener("click", (e) => {
    if (wrapper && !wrapper.contains(e.target)) {
      TOOLTIP_LOCKED = false;
      hideNodeTooltip();
    }
  });

  // klik u vieweru (prazno platno)
  if (viewer) {
    viewer.addEventListener("click", () => {
      if (TOOLTIP_LOCKED) {
        TOOLTIP_LOCKED = false;
        return; // ← NE skrivaj tooltip odmah posle klika na sferu
      }
      hideNodeTooltip();
    });
  }

  // ESC
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      hideNodeTooltip();
      closeLightbox();
    }
  });

  // lightbox
  const lightbox    = document.getElementById("lightbox");
  const lightboxImg = document.getElementById("lightboxImg");

  document.addEventListener("click", (e) => {
    if (e.target.classList.contains("zoomable")) {
      lightboxImg.src = e.target.src;
      lightbox.classList.add("open");
    }
  });

  lightbox.addEventListener("click", () => closeLightbox());
  lightboxImg.addEventListener("click", (e) => e.stopPropagation());

  function closeLightbox() {
    lightbox.classList.remove("open");
    lightboxImg.src = "";
  }
});


