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

let GRAPH_DATA = null;
let SELECTED_NODE_ID = null;
let GRAPH_LOAD_ID = 0;
let GRAPH_ABORT_CONTROLLER = null;

let LAST_ANALYSIS_KEY = null;
function analysisKey(protein, chain, type, maxdist, mindist) {
  return [protein, chain, type, maxdist, mindist].join("|");
}

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


async function loadGraph3D() {
  console.log("loadGraph3D CALLED");
  const loadId = ++GRAPH_LOAD_ID;

  // Samo poslednji klik sme da zavrsi render. Bez ovoga sporiji odgovor za
  // prethodno izabrani protein moze da stigne kasnije i pregazi novi prikaz.
  if (GRAPH_ABORT_CONTROLLER) GRAPH_ABORT_CONTROLLER.abort();
  GRAPH_ABORT_CONTROLLER = new AbortController();
  const { signal } = GRAPH_ABORT_CONTROLLER;

  const tooltip = document.getElementById("nodeTooltip");
  if (tooltip) tooltip.style.visibility = "hidden";

  const loader = document.getElementById("loader3d");
  if (loader) loader.style.visibility = "visible";

  const el = document.getElementById("viewer3d");
  if (el) el.innerHTML = "";
  GRAPH_DATA = null;
  SELECTED_NODE_ID = null;

  const proteinInput = document.getElementById("protein");
  const protein   = proteinInput?.value.trim().toUpperCase() || "";
  if (proteinInput) proteinInput.value = protein;
  const aminoname = [...document.querySelectorAll('.aa-chip.active')].map(x => x.dataset.value).join(",");
  const ss = document.getElementById("ss")?.value.trim() || "prazno";
  const sec_dist = document.getElementById("sek_dist")?.value.trim() || ""
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
    const ssDistributionMeta = document.getElementById("ssDistributionMeta");
    const ssDistributionTable = document.getElementById("ssDistributionTable");
    if (ssDistributionMeta) ssDistributionMeta.textContent = "Klikni \"Prikaži\" da učitaš podatke.";
    if (ssDistributionTable) {
      ssDistributionTable.style.display = "none";
      ssDistributionTable.querySelector("tbody").innerHTML = "";
    }
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
  loadSSDistribution(protein, chain, type, maxdist, mindist);
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
      `&ss=${encodeURIComponent(ss)}` +
      `&k=${encodeURIComponent(sec_dist)}` +
      `&type=${encodeURIComponent(type)}` +
      `&max_distance=${encodeURIComponent(maxdist)}` +
      `&min_distance=${encodeURIComponent(mindist)}`;

    console.log("Fetching graph:", url);

    const res = await fetch(url, { signal });
    if (loadId !== GRAPH_LOAD_ID) return;
    if (!res.ok) {
      throw new Error(`HTTP ${res.status} ${res.statusText}`);
    }

    const json = await res.json();
    if (loadId !== GRAPH_LOAD_ID) return;
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

    // sačuvaj sve podatke za re-rendering
    nodes.forEach(n => {
      n._x = (n.x - cx) * SCALE;
      n._y = (n.y - cy) * SCALE;
      n._z = (n.z - cz) * SCALE;
      n._aa = n.aa;
      n._index = n.index;
      n._protein = n.protein;
      n._chain = n.chain;
      n._ss = n.ss;
    });

    GRAPH_DATA = { nodes, edges, nodeById, viewer, aaColorMap3d, defaultColor };
    SELECTED_NODE_ID = null;

    const lbl = document.getElementById("viewerLabel");
    if (lbl) {
      const typeLabel = { caca: "Cα–Cα", minbezh: "Min. bez H", maxbezh: "Maks. bez H" };
      let parts = [protein.toUpperCase()];
      if (type && type !== "prazno") parts.push(typeLabel[type] || type);
      if (maxdist) parts.push(`≤ ${maxdist} Å`);
      lbl.textContent = parts.join(" · ");
      lbl.style.display = "block";
    }

    renderGraph(null);
    viewer.zoomTo();

  } catch (err) {
    if (err.name === "AbortError" || loadId !== GRAPH_LOAD_ID) return;
    console.error("loadGraph3D failed", err);
    showViewerMessage("Greška pri učitavanju grafa.", true);
  } finally {
    if (loadId === GRAPH_LOAD_ID) {
      if (loader) loader.style.visibility = "hidden";
      GRAPH_ABORT_CONTROLLER = null;
    }
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
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&include=aa_composition,sequence`;
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

async function loadSSDistribution(protein, chain, type, maxdist, mindist) {
  const meta = document.getElementById("ssDistributionMeta");
  const table = document.getElementById("ssDistributionTable");
  const tbody = table?.querySelector("tbody");

  if (!meta || !table || !tbody) return;

  meta.textContent = "Učitavanje...";
  table.style.display = "none";
  tbody.innerHTML = "";

  try {
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&include=ss_distribution`;

    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);

    const data = await res.json();
    const block = data.ss_distribution;
    const rows = block?.ss_distribution || [];
    const total = block?.total || 0;

    if (!rows.length) {
      meta.textContent = "Nema podataka o sekundarnoj strukturi za izabrani protein.";
      return;
    }

    for (const r of rows) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${r.ss ?? "–"}</td>
        <td class="num">${r.cnt ?? 0}</td>
        <td class="num">${Number(r.pct ?? 0).toFixed(1)}%</td>
      `;
      tbody.appendChild(tr);
    }

    meta.textContent = `Ukupno rezidua sa SS oznakom: ${total}`;
    table.style.display = "";
  } catch (err) {
    console.error("loadSSDistribution failed", err);
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

    const block = data.sequence;
    if (!block) return;

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

  const stamp = Date.now(); // cache-bust

  HEATMAP_URLS = {
    distance: `/api/heatmap/distance/?protein=${p}&chain=${c}&type=${t}&_=${stamp}`,
    meta: { protein, chain, type }
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
  img.alt = "Histogram dužina segmenata";
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
    if (muted) muted.textContent = "Klikni “Prikaži” da učitaš toplotnu mapu.";
    return;
  }

  // ukloni placeholder
  const muted = card.querySelector(".muted");
  if (muted) muted.remove();

  // obriši prethodne slike
  card.querySelectorAll("img.heatmap-img").forEach(x => x.remove());
  card.querySelectorAll(".heatmap-subtitle").forEach(x => x.remove());

  const { distance, meta } = HEATMAP_URLS;

  const typeLabels = {
    caca: "Cα–Cα rastojanje",
    minbezh: "Minimalno rastojanje bez H",
    maxbezh: "Maksimalno rastojanje bez H",
    prazno: "Cα–Cα rastojanje",
  };
  const typeLabel = typeLabels[meta?.type] || meta?.type || "Cα–Cα rastojanje";

  const subtitle = document.createElement("div");
  subtitle.className = "heatmap-subtitle";
  subtitle.textContent = typeLabel;
  card.appendChild(subtitle);

  const loading1 = document.createElement("div");
  loading1.className = "img-loading";
  loading1.textContent = "Učitavanje...";
  card.appendChild(loading1);

  const img1 = document.createElement("img");
  img1.className = "heatmap-img zoomable";
  img1.src = distance;
  img1.alt = `Toplotna mapa — ${typeLabel}`;
  img1.style.width = "100%";
  img1.style.borderRadius = "10px";
  img1.style.marginTop = "10px";
  img1.onload  = () => loading1.remove();
  img1.onerror = () => {
    loading1.remove();
    img1.remove();
    const err = document.createElement("div");
    err.className = "muted";
    err.textContent = "Nema podataka za toplotnu mapu (izaberi tip rastojanja i klikni Prikaži).";
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

    const fmtA = v => (v == null ? "–" : Number(v).toFixed(3) + " Å");
    const fmtP = v => (v == null ? "–" : Number(v).toFixed(1) + " %");

    const sections = [
      {
        header: "Osnovno",
        rows: [
          ["Tip",    s.type],
          ["Broj",   s.count],
        ],
      },
      {
        header: "Centralna tendencija",
        rows: [
          ["Prosek",   fmtA(s.mean)],
          ["Medijana", fmtA(s.median)],
          ["Std. devijacija", fmtA(s.std)],
        ],
      },
      {
        header: "Opseg",
        rows: [
          ["Minimum",  fmtA(s.min)],
          ["Maksimum", fmtA(s.max)],
        ],
      },
      {
        header: "Percentili",
        rows: [
          ["p10",  fmtA(s.p10)],
          ["p25",  fmtA(s.p25)],
          ["p75",  fmtA(s.p75)],
          ["p90",  fmtA(s.p90)],
          ["p95",  fmtA(s.p95)],
          ["p99",  fmtA(s.p99)],
        ],
      },
      {
        header: "Procenat kontakata ispod praga",
        rows: [
          ["< 5 Å",  fmtP(s.pct_lt_5)],
          ["< 8 Å",  fmtP(s.pct_lt_8)],
          ["< 10 Å", fmtP(s.pct_lt_10)],
          ["< 15 Å", fmtP(s.pct_lt_15)],
        ],
      },
    ];

    for (const section of sections) {
      const th = document.createElement("tr");
      th.innerHTML = `<td colspan="2" class="stat-section-header">${section.header}</td>`;
      tbody.appendChild(th);
      for (const [label, value] of section.rows) {
        const tr = document.createElement("tr");
        tr.innerHTML = `<td>${label}</td><td class="num">${value}</td>`;
        tbody.appendChild(tr);
      }
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

// =========================
//   GLOBAL STATS
// =========================
async function loadGlobalStats(attempt = 0) {
  const MAX_ATTEMPTS = 5;
  const DELAYS = [5000, 10000, 15000, 20000]; // ms između pokušaja
  const fmt = n => {
    if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + "M";
    if (n >= 1_000)     return (n / 1_000).toFixed(1) + "k";
    return String(n);
  };

  const SS_COLORS = {
    H: "#ff9f43", E: "#54a0ff", C: "#a0a0b8",
    G: "#ff6b9d", B: "#00d2d3", I: "#feca57",
    T: "#48dbfb", S: "#ff9ff3",
  };

  const setErr = (msg = "Neo4j nedostupan") => {
    ["gs-proteins","gs-aa","gs-avglen","gs-medlen","gs-minlen","gs-maxlen"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.textContent = "!";
    });
    ["gs-ss-dist","gs-top-aa","gs-dist-stats"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.textContent = msg;
    });
  };

  const retry = () => {
    if (attempt < MAX_ATTEMPTS - 1) {
      const delay = DELAYS[attempt] ?? DELAYS[DELAYS.length - 1];
      setTimeout(() => loadGlobalStats(attempt + 1), delay);
    } else {
      setErr();
    }
  };

  try {
    const res = await fetch("/api/global_stats/");
    if (!res.ok) { retry(); return; }
    const d = await res.json();
    if (d.error) { retry(); return; }

    // --- top 4 boxovi ---
    document.getElementById("gs-proteins").textContent = fmt(d.proteins    ?? 0);
    document.getElementById("gs-aa").textContent       = fmt(d.amino_acids ?? 0);
    document.getElementById("gs-avglen").textContent   = (d.avg_length    ?? "–") + " AA";
    document.getElementById("gs-medlen").textContent   = (d.median_length ?? "–") + " AA";
    document.getElementById("gs-minlen").textContent   = (d.min_length    ?? "–") + " AA";
    document.getElementById("gs-maxlen").textContent   = (d.max_length    ?? "–") + " AA";

    // --- SS distribucija ---
    const ssEl = document.getElementById("gs-ss-dist");
    if (ssEl && d.ss_distribution?.length) {
      ssEl.innerHTML = d.ss_distribution.map(s => {
        const color = SS_COLORS[s.ss] || "#888";
        return `
          <div class="gs-ss-row">
            <span class="gs-ss-label">${s.ss}</span>
            <div class="gs-ss-track">
              <div class="gs-ss-fill" style="width:${s.pct}%;background:${color}"></div>
            </div>
            <span class="gs-ss-pct">${s.pct}%</span>
          </div>`;
      }).join("");
    }

    // --- top AA ---
    const aaEl = document.getElementById("gs-top-aa");
    if (aaEl && d.top_aa?.length) {
      aaEl.innerHTML = d.top_aa.map(a =>
        `<span class="gs-aa-chip">${a.name} <b>${a.pct}%</b></span>`
      ).join("");
    }

    // --- distance stats tabela ---
    const distEl = document.getElementById("gs-dist-stats");
    if (distEl && d.dist_stats?.length) {
      const rows = d.dist_stats.map(r => `
        <tr>
          <td class="gs-dt-type">${r.type}</td>
          <td>${r.mean}</td>
          <td>${r.std}</td>
          <td>${r.min}</td>
          <td>${r.max}</td>
        </tr>`).join("");
      distEl.innerHTML = `
        <table class="gs-dist-table">
          <thead><tr>
            <th>tip</th><th>prosek</th><th>std</th><th>min</th><th>maks</th>
          </tr></thead>
          <tbody>${rows}</tbody>
        </table>`;
    }

  } catch (err) {
    console.error("loadGlobalStats failed:", err);
    retry();
  }
}

// =========================
//   GRAPH RENDERING
// =========================
function renderGraph(selectedNodeId) {
  if (!GRAPH_DATA) return;
  const { nodes, edges, nodeById, viewer, aaColorMap3d, defaultColor } = GRAPH_DATA;

  SELECTED_NODE_ID = selectedNodeId;

  // skup ID-ova direktnih suseda
  const neighborIds = new Set();
  if (selectedNodeId) {
    edges.forEach(e => {
      if (e.source === selectedNodeId) neighborIds.add(e.target);
      if (e.target === selectedNodeId) neighborIds.add(e.source);
    });
  }

  viewer.removeAllShapes();

  // --- čvorovi ---
  nodes.forEach(n => {
    const isSelected = n.id === selectedNodeId;
    const isNeighbor = neighborIds.has(n.id);
    const dimmed     = selectedNodeId && !isSelected && !isNeighbor;

    viewer.addSphere({
      center:    { x: n._x, y: n._y, z: n._z },
      radius:    isSelected ? 1.1 : 0.7,
      color:     dimmed ? 0xcccccc : (aaColorMap3d[n._aa] || defaultColor),
      opacity:   dimmed ? 0.2 : 1.0,
      clickable: true,
      callback: function(atom, v, event) {
        if (event) event.stopPropagation();
        TOOLTIP_LOCKED = true;
        const x = event?.clientX ?? NaN;
        const y = event?.clientY ?? NaN;
        showSphereTooltip(x, y, `
          <b>${n._aa}${n._index}</b><br>
          Protein: ${n._protein}<br>
          Lanac: ${n._chain}<br>
          SS: ${n._ss || "–"}
        `);
        showNodeTab(n, edges, nodeById);
        renderGraph(n.id);
      },
    });
  });

  // --- ivice ---
  edges.forEach(e => {
    const a = nodeById[e.source];
    const b = nodeById[e.target];
    if (!a || !b) return;

    const connected = selectedNodeId &&
      (e.source === selectedNodeId || e.target === selectedNodeId);
    const dimmed = selectedNodeId && !connected;

    const edgeColor = dimmed ? 0xcccccc :
      e.type === "caca"    ? 0x1f77b4 :
      e.type === "minbezh" ? 0x2ca02c :
      e.type === "maxbezh" ? 0xd62728 : 0x444444;

    viewer.addCylinder({
      start:    { x: a._x, y: a._y, z: a._z },
      end:      { x: b._x, y: b._y, z: b._z },
      radius:   connected ? 0.13 : 0.08,
      color:    edgeColor,
      opacity:  dimmed ? 0.1 : 1.0,
      clickable: !dimmed,
      callback: dimmed ? null : function(shape, v, event) {
        if (event) event.stopPropagation();
        TOOLTIP_LOCKED = true;
        const x = event?.clientX ?? NaN;
        const y = event?.clientY ?? NaN;
        showSphereTooltip(x, y, `
          <b>${a._aa}${a._index} — ${b._aa}${b._index}</b><br>
          Tip: ${e.type || "–"}<br>
          Rastojanje: ${e.value != null ? Number(e.value).toFixed(2) + " Å" : "–"}<br>
          Lanac: ${a._chain}
        `);
      },
    });
  });

  viewer.render();
}

// =========================
//   NODE TAB
// =========================
function switchToTab(tabId) {
  document.querySelectorAll(".tab-btn").forEach(b => {
    b.classList.toggle("active", b.dataset.tab === tabId);
  });
  document.querySelectorAll(".tab-content").forEach(t => {
    t.style.display = t.id === tabId ? "block" : "none";
  });
}

function showNodeTab(n, edges, nodeById) {
  // detalji čvora
  const detailCard = document.getElementById("nodeDetailCard");
  if (detailCard) {
    detailCard.innerHTML = `
      <h3>${n._aa} <span style="font-weight:400;color:var(--muted);">#${n._index}</span></h3>
      <table style="margin-top:6px;">
        <tbody>
          <tr><td>Protein</td><td class="num"><b>${n._protein}</b></td></tr>
          <tr><td>Lanac</td><td class="num"><b>${n._chain}</b></td></tr>
          <tr><td>Sek. struktura</td><td class="num"><b>${n._ss || "–"}</b></td></tr>
        </tbody>
      </table>
    `;
  }

  // susedi
  const neighbors = edges
    .filter(e => e.source === n.id || e.target === n.id)
    .map(e => {
      const neighborId = e.source === n.id ? e.target : e.source;
      const nb = nodeById[neighborId];
      return { nb, type: e.type, value: e.value };
    })
    .filter(x => x.nb)
    .sort((a, b) => (a.value ?? Infinity) - (b.value ?? Infinity));

  const neighborsCard = document.getElementById("nodeNeighborsCard");
  const countEl       = document.getElementById("nodeNeighborCount");
  const tbody         = document.querySelector("#nodeNeighborsTable tbody");
  const focusCard     = document.getElementById("nodeFocusCard");
  const focusContent  = document.getElementById("nodeFocusContent");

  if (!neighborsCard || !tbody) return;

  if (!neighbors.length) {
    neighborsCard.style.display = "none";
    if (focusCard) focusCard.style.display = "none";
  } else {
    // --- tabela kontakata ---
    countEl.textContent = `(${neighbors.length})`;
    tbody.innerHTML = "";
    for (const { nb, type, value } of neighbors) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td><b>${nb.aa}</b></td>
        <td class="num">${nb.index}</td>
        <td>${nb.chain}</td>
        <td>${nb.ss || "–"}</td>
        <td>${type || "–"}</td>
        <td class="num">${value != null ? Number(value).toFixed(3) : "–"}</td>
      `;
      tbody.appendChild(tr);
    }
    neighborsCard.style.display = "";

    // --- fokus analiza ---
    if (focusCard && focusContent) {
      const dists = neighbors.map(x => x.value).filter(v => v != null).sort((a, b) => a - b);
      let distHtml = "";
      if (dists.length) {
        const mean   = dists.reduce((s, v) => s + v, 0) / dists.length;
        const med    = dists.length % 2 === 0
          ? (dists[dists.length / 2 - 1] + dists[dists.length / 2]) / 2
          : dists[Math.floor(dists.length / 2)];
        const std    = Math.sqrt(dists.reduce((s, v) => s + (v - mean) ** 2, 0) / dists.length);
        const f = v => v.toFixed(2);
        distHtml = `
          <div class="focus-section-title">Rastojanja do suseda (Å)</div>
          <div class="focus-stat-grid">
            <div class="focus-stat-cell"><span class="fsval">${f(dists[0])}</span><span class="fslbl">Min</span></div>
            <div class="focus-stat-cell"><span class="fsval">${f(mean)}</span><span class="fslbl">Prosek</span></div>
            <div class="focus-stat-cell"><span class="fsval">${f(med)}</span><span class="fslbl">Medijana</span></div>
            <div class="focus-stat-cell"><span class="fsval">${f(std)}</span><span class="fslbl">Std</span></div>
            <div class="focus-stat-cell"><span class="fsval">${f(dists[dists.length - 1])}</span><span class="fslbl">Maks</span></div>
            <div class="focus-stat-cell"><span class="fsval">${dists.length}</span><span class="fslbl">Kontakata</span></div>
          </div>
        `;
      }

      // AA kompozicija susedstva
      const aaCounts = {};
      neighbors.forEach(({ nb }) => { aaCounts[nb.aa] = (aaCounts[nb.aa] || 0) + 1; });
      const aaSorted = Object.entries(aaCounts).sort((a, b) => b[1] - a[1]);
      const aaMax = aaSorted[0]?.[1] || 1;
      const aaHtml = `
        <div class="focus-section-title">AA kompozicija susedstva</div>
        ${aaSorted.map(([aa, cnt]) => `
          <div class="focus-bar-row">
            <span class="focus-bar-label">${aa}</span>
            <div class="focus-bar-track"><div class="focus-bar-fill" style="width:${(cnt/aaMax*100).toFixed(0)}%"></div></div>
            <span class="focus-bar-count">${cnt} (${(cnt/neighbors.length*100).toFixed(0)}%)</span>
          </div>`).join("")}
      `;

      // SS kompozicija susedstva
      const ssCounts = {};
      neighbors.forEach(({ nb }) => { const s = nb.ss || "–"; ssCounts[s] = (ssCounts[s] || 0) + 1; });
      const ssSorted = Object.entries(ssCounts).sort((a, b) => b[1] - a[1]);
      const ssMax = ssSorted[0]?.[1] || 1;
      const ssHtml = `
        <div class="focus-section-title">SS tipovi susedstva</div>
        ${ssSorted.map(([ss, cnt]) => `
          <div class="focus-bar-row">
            <span class="focus-bar-label">${ss}</span>
            <div class="focus-bar-track"><div class="focus-bar-fill ss-bar" style="width:${(cnt/ssMax*100).toFixed(0)}%"></div></div>
            <span class="focus-bar-count">${cnt} (${(cnt/neighbors.length*100).toFixed(0)}%)</span>
          </div>`).join("")}
      `;

      focusContent.innerHTML = distHtml + aaHtml + ssHtml;
      focusCard.style.display = "";
    }
  }

  switchToTab("tabNode");
}

// =========================
//   AA CHIP SELEKCIJA
// =========================
function selectAllAa() {
  document.querySelectorAll('.aa-chip').forEach(c => c.classList.add('active'));
}

function clearAa() {
  document.querySelectorAll('.aa-chip').forEach(c => c.classList.remove('active'));
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
        return;
      }
      hideNodeTooltip();
      if (SELECTED_NODE_ID !== null) {
        SELECTED_NODE_ID = null;
        renderGraph(null);
      }
    });
  }

  // ESC
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      hideNodeTooltip();
      closeLightbox();
    }
  });

  // aa chip toggle
  document.querySelectorAll('.aa-chip').forEach(chip => {
    chip.addEventListener('click', () => chip.classList.toggle('active'));
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

// Script je na kraju <body> — DOM je spreman, pozivamo odmah
loadGlobalStats();
