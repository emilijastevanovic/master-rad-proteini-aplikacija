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
let SEGMENTS_BY_SS_URL = null;

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

  // ako nemamo validne koordinate klika, tooltip ide u ugao viewera
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


// ---- dinamicka lista lanaca ------------------------------------------------
// lanci zavise od proteina, pa se lista puni tek kad se unese protein

let CHAINS_ABORT = null;
let CHAINS_LOADED_FOR = null;
let CHAINS_DEBOUNCE = null;

function setChainHint(text, cls = "") {
  const hint = document.getElementById("chainHint");
  if (!hint) return;
  hint.textContent = text;
  hint.className = "chain-hint" + (cls ? " " + cls : "");
}

function renderChainOptions(chains) {
  const sel = document.getElementById("chain");
  if (!sel) return;

  const previous = sel.value;

  sel.innerHTML = "";
  const all = document.createElement("option");
  all.value = "prazno";
  all.textContent = "– svi –";
  sel.appendChild(all);

  for (const c of chains) {
    const opt = document.createElement("option");
    opt.value = c.chain;
    opt.textContent = `${c.chain} (${c.n})`;
    sel.appendChild(opt);
  }

  // zadrzi prethodni izbor ako taj lanac postoji i u novom proteinu
  sel.value = chains.some(c => c.chain === previous) ? previous : "prazno";
}

function resetChains(hint = "", cls = "muted") {
  CHAINS_LOADED_FOR = null;
  renderChainOptions([]);
  setChainHint(hint, cls);
}

async function loadChainsFor(protein) {
  const sel = document.getElementById("chain");
  if (!sel) return;

  if (CHAINS_ABORT) CHAINS_ABORT.abort();
  CHAINS_ABORT = new AbortController();
  const { signal } = CHAINS_ABORT;

  setChainHint("…");

  try {
    const res = await fetch(`/api/chains/?protein=${encodeURIComponent(protein)}`, { signal });
    const data = await res.json();

    if (!res.ok) {
      resetChains("greška", "error");
      return;
    }

    const chains = data.chains || [];
    if (!chains.length) {
      resetChains("nema lanaca", "muted");
      return;
    }

    CHAINS_LOADED_FOR = protein;
    renderChainOptions(chains);
    setChainHint(`${chains.length}`);
  } catch (e) {
    if (e.name === "AbortError") return;
    resetChains("greška", "error");
  }
}

function refreshChains({ immediate = false } = {}) {
  const protein = document.getElementById("protein")?.value.trim().toUpperCase() || "";

  if (!protein) {
    if (CHAINS_ABORT) CHAINS_ABORT.abort();
    resetChains("");
    return;
  }

  if (protein === CHAINS_LOADED_FOR) return;

  clearTimeout(CHAINS_DEBOUNCE);
  if (immediate) {
    loadChainsFor(protein);
  } else {
    CHAINS_DEBOUNCE = setTimeout(() => loadChainsFor(protein), 350);
  }
}


async function loadGraph3D() {
  console.log("loadGraph3D CALLED");
  const loadId = ++GRAPH_LOAD_ID;

  // samo poslednji klik sme da zavrsi render, inace sporiji odgovor za
  // prethodni protein moze da pregazi novi prikaz
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

  // ako je kliknuto na Prikazi pre nego sto je debounce osvezio listu lanaca,
  // dovuci je odmah
  refreshChains({ immediate: true });

  const key = analysisKey(protein, chain, type, maxdist, mindist);
  if (key !== LAST_ANALYSIS_KEY) {
    SEQ_CACHE = null;
    SEQ_RENDERED = false;
    HEATMAP_URLS = null;
    SEGMENTS_URL = null;
    SEGMENTS_BY_SS_URL = null;
    LAST_ANALYSIS_KEY = key;

    // ocisti vidljiv sadrzaj tabova
    document.querySelectorAll(".aa-missing").forEach(x => x.remove());
    const distTable = document.getElementById("distSummaryTable");
    if (distTable) distTable.style.display = "none";
    const distMeta = document.getElementById("distSummaryMeta");
    if (distMeta) distMeta.textContent = "Klikni \"Prikaži\" da učitaš podatke.";
    const rangeTable = document.getElementById("rangeSummaryTable");
    const rangeByType = document.getElementById("rangeByTypeTable");
    const rangeByTypeTitle = document.getElementById("rangeByTypeTitle");
    const rangeMeta = document.getElementById("rangeSummaryMeta");
    const rangeFilters = document.getElementById("rangeSummaryFilters");
    if (rangeTable)  { rangeTable.style.display  = "none"; rangeTable.querySelector("tbody").innerHTML  = ""; }
    if (rangeByType) { rangeByType.style.display = "none"; rangeByType.querySelector("tbody").innerHTML = ""; }
    if (rangeByTypeTitle) rangeByTypeTitle.style.display = "none";
    if (rangeMeta) rangeMeta.textContent = "Klikni \"Prikaži\" da učitaš podatke.";
    if (rangeFilters) rangeFilters.textContent = "";
    const rangeSSPairsTable = document.getElementById("rangeSSPairsTable");
    const rangeSSPairsMeta = document.getElementById("rangeSSPairsMeta");
    const rangeSSPairsFilters = document.getElementById("rangeSSPairsFilters");
    if (rangeSSPairsTable) { rangeSSPairsTable.style.display = "none"; rangeSSPairsTable.querySelector("tbody").innerHTML = ""; }
    if (rangeSSPairsMeta) rangeSSPairsMeta.textContent = "Klikni \"Prikaži\" da učitaš podatke.";
    if (rangeSSPairsFilters) rangeSSPairsFilters.textContent = "";
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
    const closestTopMeta  = document.getElementById("closestTopMeta");
    const closestTopTable = document.getElementById("closestTopTable");
    if (closestTopMeta) closestTopMeta.textContent = "";
    if (closestTopTable) { closestTopTable.style.display = "none"; closestTopTable.querySelector("tbody").innerHTML = ""; }

    document.querySelectorAll("img.heatmap-img, img.seg-img, .img-loading").forEach(x => x.remove());
    const heatCard = document.querySelector("#tabDistance .card");
    const segCard  = document.querySelector("#tabSegments .card");
    if (heatCard && !heatCard.querySelector(".muted")) {
      const p = document.createElement("div");
      p.className = "muted";
      p.textContent = "Klikni \"Prikaži\" pa otvori ovaj tab.";
      heatCard.appendChild(p);
    }
    if (segCard) {
      segCard.querySelectorAll(".seg-subtitle, .seg-error").forEach(x => x.remove());
      if (!segCard.querySelector(".seg-placeholder")) {
        const p = document.createElement("div");
        p.className = "muted seg-placeholder";
        p.textContent = "Klikni \"Prikaži\" pa otvori ovaj tab.";
        segCard.appendChild(p);
      }
    }
  }

  // STATISTIKE
  loadAAComposition(protein, chain, type, maxdist, mindist);
  loadSSDistribution(protein, chain, type, maxdist, mindist);
  loadSSPairs(protein, chain, type, maxdist, mindist);
  loadDistanceSummary(protein, chain, type, maxdist, mindist);
  loadRangeSummary(protein, chain, type, maxdist, mindist, ss, sec_dist, aminoname);
  loadOutlierSS(protein, chain, type, maxdist, mindist);

  prefetchSeq(protein, chain, type, maxdist, mindist);

  setHeatmapUrls(protein, chain, type, maxdist, mindist, aminoname);
  setSegmentsUrl(protein, chain);

  // ako je tab vec otvoren, odmah re-renderuj
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
    // encodeURIComponent za svaki parametar
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
      // 400 nosi objasnjenje sta je pogresno uneto, prikazi ga korisniku
      let detail = "";
      try { detail = (await res.json())?.detail || ""; } catch (_) {}
      const httpErr = new Error(detail || `HTTP ${res.status} ${res.statusText}`);
      httpErr.userMessage = detail;
      throw httpErr;
    }

    const json = await res.json();
    if (loadId !== GRAPH_LOAD_ID) return;
    let nodes = json.nodes || [];
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

    // cvor bez sve tri koordinate se ne sme pustiti dalje: (null - cx) daje
    // NaN, a jedan NaN u zoomTo() obori ceo prikaz
    const imaKoordinate = n =>
      Number.isFinite(n.x) && Number.isFinite(n.y) && Number.isFinite(n.z);
    const bezKoordinata = nodes.length - nodes.filter(imaKoordinate).length;
    nodes = nodes.filter(imaKoordinate);

    if (!nodes.length) {
      showViewerMessage(
        "Nijedan reziduum za izabrane filtere nema Cα koordinate, pa 3D prikaz nije moguć.",
        true
      );
      viewer.render();
      return;
    }
    if (bezKoordinata) {
      console.warn(`${bezKoordinata} čvorova bez koordinata — izostavljeni iz 3D prikaza.`);
    }

    const xs = nodes.map(n => n.x);
    const ys = nodes.map(n => n.y);
    const zs = nodes.map(n => n.z);

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
      // nepolarne (hidrofobne), zuta
      ALA: 0xe6b800, VAL: 0xe6b800, ILE: 0xe6b800, LEU: 0xe6b800,
      MET: 0xe6b800, PHE: 0xe6b800, TRP: 0xe6b800, PRO: 0xe6b800,
      GLY: 0xe6b800,
      // polarne, zelena
      SER: 0x4caf50, THR: 0x4caf50, CYS: 0x4caf50, TYR: 0x4caf50,
      ASN: 0x4caf50, GLN: 0x4caf50,
      // pozitivne, plava
      LYS: 0x2979ff, ARG: 0x2979ff, HIS: 0x2979ff,
      // negativne, crvena
      ASP: 0xe53935, GLU: 0xe53935,
    };
    const defaultColor = 0x9e9e9e;

    // sacuvaj sve podatke za re-rendering
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
      if (mindist && maxdist) parts.push(`${mindist}–${maxdist} Å`);
      else if (mindist) parts.push(`≥ ${mindist} Å`);
      else if (maxdist) parts.push(`≤ ${maxdist} Å`);
      lbl.textContent = parts.join(" · ");
      lbl.style.display = "block";
    }

    renderGraph(null);
    viewer.zoomTo();

  } catch (err) {
    if (err.name === "AbortError" || loadId !== GRAPH_LOAD_ID) return;
    console.error("loadGraph3D failed", err);
    showViewerMessage(err.userMessage || "Greška pri učitavanju grafa.", true);
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

    // ocekujemo: data.aa_composition = { protein, total, rows }
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

      // procenat je vec izracunat u Cypheru (round(cnt*100/total, 1))
      const pct = Number(r.pct ?? 0);

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

  // ako jos nista nije prefetched
  if (!SEQ_CACHE) {
    seqEl.textContent = "Klikni “Prikaži” da učitaš podatke.";
    if (ssEl) ssEl.textContent = "";
    return;
  }

  // ako je vec renderovano, samo izadji
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

function setHeatmapUrls(protein, chain, type, maxdist, mindist, aminoname) {
  const p = encodeURIComponent(protein);
  const c = encodeURIComponent(chain);
  const t = encodeURIComponent(type);

  const stamp = Date.now(); // cache-bust

  // druga mapa se crta samo za pravi podskup aminokiselina, inace bi bila ista kao prva
  const selected = (aminoname || "").split(",").filter(Boolean);
  const total = document.querySelectorAll(".aa-chip").length;
  const isSubset = selected.length > 0 && selected.length < total;

  // treca mapa je agregirana po vrsti aminokiseline, pa ima smisla tek od dve navise
  const aaSelekcija = encodeURIComponent(selected.join(","));

  HEATMAP_URLS = {
    distance: `/api/heatmap/distance/?protein=${p}&chain=${c}&type=${t}&_=${stamp}`,
    subset: isSubset
      ? `/api/heatmap/distance/?protein=${p}&chain=${c}&type=${t}` +
        `&aminoname=${aaSelekcija}&_=${stamp}`
      : null,
    aaTypes: selected.length >= 2
      ? `/api/heatmap/aa_types/?protein=${p}&chain=${c}&type=${t}` +
        `&aminoname=${aaSelekcija}&_=${stamp}`
      : null,
    meta: { protein, chain, type, aminonames: selected, isSubset }
  };
}

function setSegmentsUrl(protein, chain) {
  const stamp = Date.now();
  const q = `protein=${encodeURIComponent(protein)}&chain=${encodeURIComponent(chain)}&_=${stamp}`;
  SEGMENTS_URL = `/api/heatmap/segments/?${q}`;
  SEGMENTS_BY_SS_URL = `/api/heatmap/segments_by_ss/?${q}`;
}

function renderSegmentsImg() {
  const tab = document.getElementById("tabSegments");
  if (!tab) return;
  const card = tab.querySelector(".card");
  if (!card) return;

  if (!SEGMENTS_URL) {
    const placeholder = card.querySelector(".seg-placeholder");
    if (placeholder) placeholder.textContent = "Klikni \"Prikaži\" da učitaš histogram.";
    return;
  }

  // brise se samo poruka za klik i greske iz prethodnog pokusaja,
  // objasnjenje sta je segment ostaje u kartici
  card.querySelectorAll(".seg-placeholder, .seg-error").forEach(x => x.remove());
  card.querySelectorAll("img.seg-img").forEach(x => x.remove());
  card.querySelectorAll(".seg-subtitle").forEach(x => x.remove());

  const addImg = (url, alt, errText) => {
    const loading = document.createElement("div");
    loading.className = "img-loading";
    loading.textContent = "Učitavanje...";
    card.appendChild(loading);

    const img = document.createElement("img");
    img.className = "seg-img zoomable";
    img.src = url;
    img.alt = alt;
    img.style.width = "100%";
    img.style.borderRadius = "10px";
    img.style.marginTop = "10px";
    img.onload  = () => loading.remove();
    img.onerror = () => {
      loading.remove();
      img.remove();
      const err = document.createElement("div");
      err.className = "muted seg-error";
      err.textContent = errText;
      card.appendChild(err);
    };
    card.appendChild(img);
  };

  const sub1 = document.createElement("div");
  sub1.className = "heatmap-subtitle seg-subtitle";
  sub1.textContent = "Svi SS tipovi";
  card.appendChild(sub1);

  addImg(SEGMENTS_URL, "Histogram dužina segmenata", "Nema podataka za segmente (klikni Prikaži).");

  if (!SEGMENTS_BY_SS_URL) return;

  const sub2 = document.createElement("div");
  sub2.className = "heatmap-subtitle seg-subtitle";
  sub2.style.marginTop = "18px";
  sub2.textContent = "Po SS tipu";
  card.appendChild(sub2);

  addImg(
    SEGMENTS_BY_SS_URL,
    "Histogrami dužina segmenata po SS tipu",
    "Nema podataka o sekundarnoj strukturi za ovaj protein/lanac."
  );
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

  // obrisi prethodne slike
  card.querySelectorAll("img.heatmap-img").forEach(x => x.remove());
  card.querySelectorAll(".heatmap-subtitle").forEach(x => x.remove());

  card.querySelectorAll(".heatmap-err").forEach(x => x.remove());

  const { distance, subset, aaTypes, meta } = HEATMAP_URLS;

  const typeLabels = {
    caca: "Cα–Cα rastojanje",
    minbezh: "Minimalno rastojanje bez H",
    maxbezh: "Maksimalno rastojanje bez H",
    prazno: "Cα–Cα rastojanje",
  };
  const typeLabel = typeLabels[meta?.type] || meta?.type || "Cα–Cα rastojanje";
  const aaList = (meta?.aminonames || []).join(", ");

  // tri mape se razlikuju samo po izvoru, naslovu i poruci o gresci
  const dodajMapu = ({ src, naslov, alt, greska, prva }) => {
    const subtitle = document.createElement("div");
    subtitle.className = "heatmap-subtitle";
    if (!prva) subtitle.style.marginTop = "18px";
    subtitle.textContent = naslov;
    card.appendChild(subtitle);

    const loading = document.createElement("div");
    loading.className = "img-loading";
    loading.textContent = "Učitavanje...";
    card.appendChild(loading);

    const img = document.createElement("img");
    img.className = "heatmap-img zoomable";
    img.src = src;
    img.alt = alt;
    img.style.width = "100%";
    img.style.borderRadius = "10px";
    img.style.marginTop = "10px";
    img.onload  = () => loading.remove();
    img.onerror = () => {
      loading.remove();
      img.remove();
      const err = document.createElement("div");
      err.className = "muted heatmap-err";
      err.textContent = greska;
      card.appendChild(err);
    };
    card.appendChild(img);
  };

  dodajMapu({
    src: distance,
    naslov: typeLabel,
    alt: `Toplotna mapa — ${typeLabel}`,
    greska: "Nema podataka za toplotnu mapu (izaberi tip rastojanja i klikni Prikaži).",
    prva: true,
  });

  // druga mapa: samo oznacene aminokiseline, svaki reziduum sa svakim
  if (subset) {
    dodajMapu({
      src: subset,
      naslov: `${typeLabel}: selektovane ${aaList}`,
      alt: `Toplotna mapa — ${aaList}`,
      greska: `Nema rastojanja između selektovanih aminokiselina (${aaList}).`,
    });
  }

  // treca mapa: prosek po vrsti, matrica je velika koliko ima oznacenih vrsta
  if (aaTypes) {
    dodajMapu({
      src: aaTypes,
      naslov: `${typeLabel}: prosek po vrsti aminokiseline (${aaList})`,
      alt: `Prosečno rastojanje po vrsti aminokiseline — ${aaList}`,
      greska: `Nema parova za prosek po vrsti (${aaList}).`,
    });
  }
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

const DIST_TYPE_LABELS = {
  caca: "Cα–Cα",
  minbezh: "Min. bez H",
  maxbezh: "Maks. bez H",
};

// opisni red iznad tabela: tip rastojanja i izabrani opseg
// lanac, k, ss i aminokiseline se dopisuju samo ako su zadati, da ostane jedan red
function describeRangeFilters(f) {
  if (!f) return "";

  const type = f.type ? (DIST_TYPE_LABELS[f.type] || f.type) : "Sva";

  let range;
  if (f.min_distance != null && f.max_distance != null) {
    range = `u opsegu ${f.min_distance}–${f.max_distance} Å`;
  } else if (f.max_distance != null) {
    range = `do ${f.max_distance} Å`;
  } else if (f.min_distance != null) {
    range = `od ${f.min_distance} Å`;
  } else {
    range = "bez zadatih granica";
  }

  const extra = [];
  if (f.chain) extra.push(`lanac ${f.chain}`);
  if (f.k != null) extra.push(`k ≥ ${f.k}`);
  if (f.ss) extra.push(`SS ${f.ss}`);
  if (f.aminonames && f.aminonames.length && f.aminonames.length < 20) {
    extra.push(`${f.aminonames.length} od 20 AA`);
  }

  return `${type} rastojanja ${range}` + (extra.length ? ` (${extra.join(", ")})` : "") + ".";
}

// tabela ss parova u izabranom opsegu, puni se iz istog odgovora kao i pregled opsega
function renderRangeSSPairs(s) {
  const meta      = document.getElementById("rangeSSPairsMeta");
  const filtersEl = document.getElementById("rangeSSPairsFilters");
  const table     = document.getElementById("rangeSSPairsTable");
  const tbody     = table?.querySelector("tbody");

  if (!meta || !table || !tbody) return;

  tbody.innerHTML = "";
  table.style.display = "none";

  if (!s) {
    meta.textContent = "Nema podataka.";
    if (filtersEl) filtersEl.textContent = "";
    return;
  }

  if (filtersEl) filtersEl.textContent = describeRangeFilters(s.filters);

  const rows = s.by_ss_pair || [];
  if (!rows.length) {
    meta.textContent = "Nema parova sa poznatom sekundarnom strukturom u izabranom opsegu.";
    return;
  }

  const fmt3 = v => (v == null ? "–" : Number(v).toFixed(3));

  for (const r of rows) {
    const tr = document.createElement("tr");
    tr.innerHTML =
      `<td>${r.ss_pair}</td>` +
      `<td class="num">${Number(r.count).toLocaleString("sr-RS")}</td>` +
      `<td class="num">${r.pct == null ? "–" : Number(r.pct).toFixed(1)}</td>` +
      `<td class="num">${fmt3(r.mean)}</td>` +
      `<td class="num">${fmt3(r.median)}</td>` +
      `<td class="num">${fmt3(r.min)}</td>` +
      `<td class="num">${fmt3(r.max)}</td>` +
      `<td class="num">${fmt3(r.std)}</td>`;
    tbody.appendChild(tr);
  }

  // filter po ss trazi da obe aminokiseline budu iste strukture, pa tabela
  // tada ima samo jedan red
  meta.textContent = s.filters?.ss
    ? `Izabran je filter SS „${s.filters.ss}", pa ulaze samo parovi kod kojih obe aminokiseline imaju tu strukturu.`
    : "";

  table.style.display = "";
}

async function loadRangeSummary(protein, chain, type, maxdist, mindist, ss, k, aminoname) {
  const meta       = document.getElementById("rangeSummaryMeta");
  const filtersEl  = document.getElementById("rangeSummaryFilters");
  const table      = document.getElementById("rangeSummaryTable");
  const tbody      = table?.querySelector("tbody");
  const typeTitle  = document.getElementById("rangeByTypeTitle");
  const typeTable  = document.getElementById("rangeByTypeTable");
  const typeTbody  = typeTable?.querySelector("tbody");
  const ssPairsMeta = document.getElementById("rangeSSPairsMeta");

  if (!meta || !table || !tbody || !typeTable || !typeTbody) return;

  meta.textContent = "Učitavanje...";
  if (filtersEl) filtersEl.textContent = "";
  if (ssPairsMeta) ssPairsMeta.textContent = "Učitavanje...";
  table.style.display = "none";
  tbody.innerHTML = "";
  typeTable.style.display = "none";
  if (typeTitle) typeTitle.style.display = "none";
  typeTbody.innerHTML = "";

  try {
    const url =
      `/api/stats/?protein=${encodeURIComponent(protein)}` +
      `&chain=${encodeURIComponent(chain)}` +
      `&type=${encodeURIComponent(type)}` +
      `&max_distance=${encodeURIComponent(maxdist)}` +
      `&min_distance=${encodeURIComponent(mindist)}` +
      `&ss=${encodeURIComponent(ss || "")}` +
      `&k=${encodeURIComponent(k || "")}` +
      `&aminoname=${encodeURIComponent(aminoname || "")}` +
      `&include=range_summary`;

    const res = await fetch(url);
    if (!res.ok) {
      const detail = await res.json().catch(() => null);
      const msg = detail?.detail || `Greška pri učitavanju (HTTP ${res.status}).`;
      meta.textContent = msg;
      if (ssPairsMeta) ssPairsMeta.textContent = msg;
      return;
    }

    const s = (await res.json()).range_summary;
    renderRangeSSPairs(s);

    if (!s) {
      meta.textContent = "Nema podataka.";
      return;
    }

    if (filtersEl) filtersEl.textContent = describeRangeFilters(s.filters);

    if (s.empty) {
      meta.textContent = s.total_base
        ? `Nijedan od ${s.total_base.toLocaleString("sr-RS")} parova ne upada u izabrani opseg.`
        : "Nema parova za izabrane filtere.";
      return;
    }

    const fmtA = v => (v == null ? "–" : Number(v).toFixed(3) + " Å");
    const fmtN = v => (v == null ? "–" : Number(v).toLocaleString("sr-RS"));

    const sections = [
      {
        header: "Broj parova",
        rows: [
          ["U izabranom opsegu", fmtN(s.count)],
          ["Bez granica", fmtN(s.total_base)],
          ["Udeo", s.pct_of_base == null ? "–" : Number(s.pct_of_base).toFixed(1) + " %"],
        ],
      },
      {
        header: "Centralna tendencija",
        rows: [
          ["Prosek", fmtA(s.mean)],
          ["Medijana", fmtA(s.median)],
          ["Std. devijacija", fmtA(s.std)],
        ],
      },
      {
        header: "Opseg vrednosti",
        rows: [
          ["Minimum", fmtA(s.min)],
          ["p25", fmtA(s.p25)],
          ["p75", fmtA(s.p75)],
          ["Maksimum", fmtA(s.max)],
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

    const byType = s.by_type || [];
    if (byType.length) {
      for (const r of byType) {
        const tr = document.createElement("tr");
        tr.innerHTML =
          `<td>${DIST_TYPE_LABELS[r.type] || r.type}</td>` +
          `<td class="num">${fmtN(r.count)}</td>` +
          `<td class="num">${r.mean == null ? "–" : Number(r.mean).toFixed(3)}</td>` +
          `<td class="num">${r.median == null ? "–" : Number(r.median).toFixed(3)}</td>` +
          `<td class="num">${r.min == null ? "–" : Number(r.min).toFixed(3)}</td>` +
          `<td class="num">${r.max == null ? "–" : Number(r.max).toFixed(3)}</td>` +
          `<td class="num">${r.std == null ? "–" : Number(r.std).toFixed(3)}</td>`;
        typeTbody.appendChild(tr);
      }
      typeTable.style.display = "";
      if (typeTitle) typeTitle.style.display = "";
    }

  } catch (err) {
    console.error("loadRangeSummary failed", err);
    meta.textContent = "Greška pri učitavanju (proveri konzolu).";
    if (ssPairsMeta) ssPairsMeta.textContent = "Greška pri učitavanju (proveri konzolu).";
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
  const closeMeta  = document.getElementById("closestTopMeta");
  const closeTable = document.getElementById("closestTopTable");
  const closeTbody = closeTable?.querySelector("tbody");

  if (!ssMeta || !ssTable || !ssTbody || !topMeta || !topTable || !topTbody) return;

  ssMeta.textContent  = "Učitavanje...";
  topMeta.textContent = "";
  ssTable.style.display  = "none";
  topTable.style.display = "none";
  ssTbody.innerHTML  = "";
  topTbody.innerHTML = "";
  if (closeMeta)  closeMeta.textContent = "";
  if (closeTable) closeTable.style.display = "none";
  if (closeTbody) closeTbody.innerHTML = "";

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
    const closestResidues = block?.closest_residues || [];

    const fmt = v => (v == null ? "–" : Number(v).toFixed(3));

    const residueRow = r => `
      <td>${r.aa ?? "–"}</td>
      <td class="num">${r.index ?? "–"}</td>
      <td>${r.chain ?? "–"}</td>
      <td>${r.ss ?? "–"}</td>
      <td class="num">${fmt(r.avg_dist)}</td>
    `;

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

    if (closestResidues.length && closeTbody && closeTable && closeMeta) {
      for (const r of closestResidues) {
        const tr = document.createElement("tr");
        tr.innerHTML = residueRow(r);
        closeTbody.appendChild(tr);
      }
      // meta se popunjava samo kad protein ima manje od 15 aminokiselina
      closeMeta.textContent = closestResidues.length < 15
        ? `Protein ima svega ${closestResidues.length} aminokiselina.` : "";
      closeTable.style.display = "";
    }

    if (topResidues.length) {
      for (const r of topResidues) {
        const tr = document.createElement("tr");
        tr.innerHTML = residueRow(r);
        topTbody.appendChild(tr);
      }
      topMeta.textContent = "";
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
  const DELAYS = [5000, 10000, 15000, 20000]; // ms izmedju pokusaja
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

  // --- cvorovi ---
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
  // detalji cvora
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
  if (ss === "G") return "ss-G";
  if (ss === "I") return "ss-I";
  if (ss === "E") return "ss-E";
  if (ss === "B") return "ss-B";
  if (ss === "T") return "ss-T";
  if (ss === "S") return "ss-S";
  return "ss-dot";
}

// DSSP oznake, redosled je i redosled u legendi
const SS_LEGEND = [
  ["H", "α-heliks"],
  ["G", "3₁₀-heliks"],
  ["I", "π-heliks"],
  ["E", "β-lanac"],
  ["B", "β-most"],
  ["T", "okret"],
  ["S", "savijanje"],
  ["-", "nesvrstano"],
];

function renderSsLegend() {
  const box = document.getElementById("ssLegend");
  if (!box || box.childElementCount) return;

  for (const [code, label] of SS_LEGEND) {
    const item = document.createElement("div");
    item.className = "ss-legend-item";
    item.innerHTML =
      `<span class="ss-legend-key ${ssToClass(code)}">${code}</span><span>${label}</span>`;
    box.appendChild(item);
  }
}


const AA_HYDROPHOBIC = new Set(["ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO","GLY"]);
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

// samo znacenje boja, koja je aminokiselina u kojoj grupi vidi se iz obojenog niza
const AA_LEGEND = [
  ["k-nepolarne", "nepolarne (hidrofobne)"],
  ["k-polarne",   "polarne"],
  ["k-bazne",     "bazne (pozitivne)"],
  ["k-kisele",    "kisele (negativne)"],
  ["k-ostalo",    "ostalo"],
];

function renderAaLegend() {
  const box = document.getElementById("aaLegend");
  if (!box || box.childElementCount) return;

  for (const [swatch, label] of AA_LEGEND) {
    const item = document.createElement("div");
    item.className = "aa-legend-item";
    item.innerHTML = `<span class="aa-legend-key ${swatch}"></span><span>${label}</span>`;
    box.appendChild(item);
  }
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

      // 3) prikazi izabrani
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

// =========================
//   ŠIRINA PANELA SA STATISTIKAMA
// =========================
const STATS_WIDTH_KEY = "statsPanelWidth";
const STATS_WIDTH_DEFAULT = 380;
const STATS_WIDTH_MIN = 300;
const GRAPH_WIDTH_MIN = 260;   // koliko najmanje ostaje 3D prikazu

// 3Dmol platno ne prati promenu sirine samo od sebe
function resizeViewer() {
  const viewer = GRAPH_DATA && GRAPH_DATA.viewer;
  if (!viewer || typeof viewer.resize !== "function") return;
  requestAnimationFrame(() => {
    viewer.resize();
    viewer.render();
  });
}

function setStatsWidth(layout, px) {
  const maxW = Math.max(STATS_WIDTH_MIN, layout.clientWidth - GRAPH_WIDTH_MIN);
  const w = Math.round(Math.min(maxW, Math.max(STATS_WIDTH_MIN, px)));
  layout.style.setProperty("--stats-w", w + "px");
  return w;
}

function setupStatsResize() {
  const layout   = document.querySelector(".layout");
  const splitter = document.getElementById("statsSplitter");
  const fullBtn  = document.getElementById("statsFullBtn");
  if (!layout) return;

  // zapamcena sirina iz prethodne sesije
  const saved = parseInt(localStorage.getItem(STATS_WIDTH_KEY) || "", 10);
  if (Number.isFinite(saved)) setStatsWidth(layout, saved);

  if (splitter) {
    splitter.addEventListener("pointerdown", (e) => {
      e.preventDefault();
      splitter.setPointerCapture(e.pointerId);
      document.body.classList.add("is-resizing-stats");

      const onMove = (ev) => {
        // panel je desno: sirina = desna ivica layout-a minus pozicija kursora
        const rect = layout.getBoundingClientRect();
        setStatsWidth(layout, rect.right - ev.clientX);
      };

      const onUp = () => {
        splitter.removeEventListener("pointermove", onMove);
        splitter.removeEventListener("pointerup", onUp);
        splitter.removeEventListener("pointercancel", onUp);
        document.body.classList.remove("is-resizing-stats");
        const cur = layout.style.getPropertyValue("--stats-w");
        if (cur) localStorage.setItem(STATS_WIDTH_KEY, parseInt(cur, 10));
        resizeViewer();
      };

      splitter.addEventListener("pointermove", onMove);
      splitter.addEventListener("pointerup", onUp);
      splitter.addEventListener("pointercancel", onUp);
    });

    // dupli klik vraca podrazumevanu sirinu
    splitter.addEventListener("dblclick", () => {
      setStatsWidth(layout, STATS_WIDTH_DEFAULT);
      localStorage.removeItem(STATS_WIDTH_KEY);
      resizeViewer();
    });
  }

  const graphFullBtn = document.getElementById("graphFullBtn");

  // cela sirina za panel i za graf, iskljucuju jedno drugo
  function applyFullMode(mode) {   // mode: "stats" | "graph" | null
    layout.classList.toggle("stats-full", mode === "stats");
    layout.classList.toggle("graph-full", mode === "graph");

    if (fullBtn) {
      fullBtn.textContent = mode === "stats" ? "⤡ Vrati prikaz" : "⤢ Cela širina";
      fullBtn.title = mode === "stats"
        ? "Vrati 3D prikaz i uobičajenu širinu panela"
        : "Panel preko cele širine — pogodno za screenshot tabela";
    }
    if (graphFullBtn) {
      graphFullBtn.textContent = mode === "graph" ? "⤡ Vrati panel" : "⤢ Cela širina";
      graphFullBtn.title = mode === "graph"
        ? "Vrati panel sa statistikama"
        : "Graf preko cele širine — pogodno za screenshot 3D prikaza";
    }

    resizeViewer();
  }

  if (fullBtn) {
    fullBtn.addEventListener("click", () => {
      applyFullMode(layout.classList.contains("stats-full") ? null : "stats");
    });
  }

  if (graphFullBtn) {
    graphFullBtn.addEventListener("click", () => {
      applyFullMode(layout.classList.contains("graph-full") ? null : "graph");
    });
  }

  // kod suzavanja prozora ne dozvoli da panel proguta ceo prostor
  window.addEventListener("resize", () => {
    const cur = parseInt(layout.style.getPropertyValue("--stats-w"), 10);
    if (Number.isFinite(cur)) setStatsWidth(layout, cur);
  });
}

function hideNodeTooltip() {
  const tooltip = document.getElementById("nodeTooltip");
  if (tooltip) tooltip.style.visibility = "hidden";
}

// setupTabs se poziva kad se DOM ucita
document.addEventListener("DOMContentLoaded", () => {
  // 1) inicijalizuj tabove
  setupTabs();
  setupStatsResize();
  renderSsLegend();
  renderAaLegend();

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

  // lanci se ucitavaju cim se unese protein
  const proteinInput = document.getElementById("protein");
  if (proteinInput) {
    proteinInput.addEventListener("input", () => refreshChains());
    proteinInput.addEventListener("change", () => refreshChains({ immediate: true }));
    refreshChains({ immediate: true });   // ako je polje vec popunjeno (npr. posle reload-a)
  }

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

// script je na kraju <body>, DOM je spreman pa zovemo odmah
loadGlobalStats();
