let SEQ_CACHE = null;     
let SEQ_RENDERED = false;

let HEATMAP_URLS = null;

let LAST_ANALYSIS_KEY = null;
function analysisKey(protein, chain, type, maxdist, mindist) {
  return [protein, chain, type, maxdist, mindist].join("|");
}


async function loadGraph() {
  const protein = document.getElementById("protein").value.trim();
  const type = document.getElementById("type").value;
  const maxdist = document.getElementById("maxdist").value.trim();

  if (!protein) {
    alert("Unesi naziv proteina!");
    return;
  }

  const url = `/api/graph/?protein=${protein}&type=${type}&max_distance=${maxdist}`;
  const response = await fetch(url);
  const data = await response.json();

  console.log("Podaci iz API-ja:", data);

  // ✅ API sada vraća {nodes: [...], edges: [...]}
  const nodesData = data.nodes || [];
  const edgesData = data.edges || [];

  // Pretvaranje u Cytoscape elemente
  const nodes = nodesData.map(n => ({
    data: {
      id: n.id,
      label: `${n.aa}${n.index}`,
    },
  }));

  const edges = edgesData.map(e => ({
    data: {
      id: `${e.source}-${e.target}`,
      source: e.source,
      target: e.target,
      label: e.value ? e.value.toFixed(2) : "",
      type: e.type,
    },
  }));

  const elements = [...nodes, ...edges];

  //  Ako postoji prethodni graf — očisti ga pre novog učitavanja
  const container = document.getElementById("graph");
  container.innerHTML = "";

  // Prikaz grafa
  const cy = cytoscape({
    container: container,
    elements: elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "#0074D9",
          "label": "data(label)",
          "color": "#fff",
          "text-valign": "center",
          "text-outline-width": 1,
          "text-outline-color": "#0074D9",
          "font-size": 10,
        },
      },
      {
        selector: "edge",
        style: {
          "width": 1.5,
          "line-color": "#aaa",
          "curve-style": "bezier",
          "label": "data(label)",
          "font-size": 8,
          "color": "#555",
        },
      },
      { selector: 'edge[type="caca"]', style: { "line-color": "#1f77b4" } },
      { selector: 'edge[type="minbezh"]', style: { "line-color": "#2ca02c" } },
      { selector: 'edge[type="maxbezh"]', style: { "line-color": "#d62728" } },
    ],
    layout: {
      name: "cose",
      padding: 20,
    },
    wheelSensitivity: 0.2,
  });
}


// =========================
//   SHOW TOOLTIP FUNCTION
// =========================


function showSphereTooltip(pixelX, pixelY, html) {
  const tooltip = document.getElementById("nodeTooltip");

  tooltip.innerHTML = html;
  tooltip.style.left = (pixelX + 12) + "px";
  tooltip.style.top = (pixelY + 12) + "px";
  tooltip.style.visibility = "visible";
}

// =========================
//   3D GRAPH
// =========================
async function loadGraph3D() {
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
    alert("Unesi naziv proteina!");
    if (loader) loader.style.visibility = "hidden";
    return;
  }

  const key = analysisKey(protein, chain, type, maxdist, mindist);
  if (key !== LAST_ANALYSIS_KEY) {
    // reset caches za novu analizu
    SEQ_CACHE = null;
    SEQ_RENDERED = false;
    HEATMAP_URLS = null;
    LAST_ANALYSIS_KEY = key;
  }

  // STATISTIKE
  if (typeof loadAAComposition === "function") {
    loadAAComposition(protein, chain, type, maxdist, mindist);
  }
  
  prefetchSeq(protein, chain, type, maxdist, mindist);
  
  // ✅ samo formiramo URL-e za img src
  setHeatmapUrls(protein, chain, type, maxdist, mindist);

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

    const ssColorMap = {
      H: 0xff7f0e,
      E: 0x2ca02c,
      T: 0xffd700,
      S: 0x9467bd
    };
    const defaultColor = 0x1f77b4;

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
        color: ssColorMap[n.ss] || defaultColor,
        clickable: true,
        callback: function(atom, viewer, event) {
          if (event) event.stopPropagation();
          const x = event?.clientX ?? 0;
          const y = event?.clientY ?? 0;

          const html = `
            <b>${n._aa}${n._index}</b><br>
            Protein: ${n._protein}<br>
            Lanac: ${n._chain}<br>
            SS: ${n._ss || "-"}
          `;
          if (typeof showSphereTooltip === "function") {
            showSphereTooltip(x, y, html);
          }
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

      viewer.addLine({
        start: { x: a._x, y: a._y, z: a._z },
        end:   { x: b._x, y: b._y, z: b._z },
        color: edgeColor,
        linewidth: 1.0
      });
    });

    viewer.zoomTo();
    viewer.render();

  } catch (err) {
    console.error("loadGraph3D failed ❌", err);
    alert("Greška pri učitavanju grafa. Pogledaj Console.");
  } finally {
    if (loader) loader.style.visibility = "hidden";
  }
}

async function loadGraph3D1() {

  const tooltip = document.getElementById("nodeTooltip");
  tooltip.style.visibility = "hidden";

  const loader = document.getElementById("loader3d");
  loader.style.visibility = "visible";

  const el = document.getElementById("viewer3d");
  el.innerHTML = "";

  const protein   = document.getElementById("protein").value.trim();
  const aminoname = document.getElementById("aminoname").value.trim();
  const type      = document.getElementById("type").value;
  const maxdist   = document.getElementById("maxdist").value.trim();
  const mindist   = document.getElementById("mindist").value.trim();
  const chain     = document.getElementById("chain").value.trim();

  if (!protein) {
      alert("Unesi naziv proteina!");
      return;
  }

  const url = `/api/graph3d/?protein=${protein}&aminoname=${aminoname}&chain=${chain}&type=${type}&max_distance=${maxdist}&min_distance=${mindist}`;
  const res = await fetch(url);
  const { nodes, edges } = await res.json();

  const viewer = $3Dmol.createViewer(el, { backgroundColor: "white" });

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

  const ssColorMap = {
      H: 0xff7f0e,
      E: 0x2ca02c,
      T: 0xffd700,
      S: 0x9467bd
  };

  const defaultColor = 0x1f77b4;

  // --------------------
  // Draw spheres
  // --------------------
  nodes.forEach(n => {
      n._x = (n.x - cx)*SCALE;
      n._y = (n.y - cy)*SCALE;
      n._z = (n.z - cz)*SCALE;

      n._aa = n.aa
      n._index = n.index
      n._protein = n.protein
      n._chain = n.chain
      n._ss = n.ss

      viewer.addSphere({
          center: { x: n._x, y: n._y, z: n._z },
          radius: 0.7,
          color: ssColorMap[n.ss] || defaultColor,
          clickable: true,


          callback: function(atom, viewer, event) {

              event.stopPropagation();
              const x = event.clientX;
              const y = event.clientY;

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

  // --------------------
  // Draw edges
  // --------------------
  edges.forEach(e => {
      const a = nodeById[e.source];
      const b = nodeById[e.target];
      if (!a || !b) return;

      const edgeColor =
          e.type === "caca"    ? 0x1f77b4 :
          e.type === "minbezh" ? 0x2ca02c :
          e.type === "maxbezh" ? 0xd62728 :
                                 0x444444;

      viewer.addLine({
          start: { x: a._x, y: a._y, z: a._z },
          end:   { x: b._x, y: b._y, z: b._z },
          color: edgeColor,
          linewidth: 1.0
      });
  });


  viewer.zoomTo();
  viewer.render();
  loader.style.visibility = "hidden";

  if (protein) {
    loadAAComposition(protein, chain, type, maxdist, mindist);
  }
}

document.addEventListener("click", function (e) {
  const tooltip = document.getElementById("nodeTooltip");
  const wrapper = document.getElementById("viewer3d-wrapper");

  if (!tooltip || !wrapper) return;

  if (!wrapper.contains(e.target)) {
    tooltip.style.visibility = "hidden";
  }
});



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

    meta.textContent = `Ukupno aminokiselina: ${total} · Različitih tipova: ${rows.length}`;
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
    span.className = `seqToken ${ssToClass(ss)}`;
    span.textContent = `[${aa}]`;

    // tooltip
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
    distance: `/api/heatmap/distance?protein=${p}&chain=${c}&type=${t}&max_distance=${maxd}&min_distance=${mind}&_=${stamp}`,
    meta: { protein, chain, type, maxdist, mindist }
  };
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

  // distance img
  const img1 = document.createElement("img");
  img1.className = "heatmap-img";
  img1.src = distance;
  img1.alt = "Distance heatmap";
  img1.loading = "lazy";
  img1.style.width = "100%";
  img1.style.borderRadius = "10px";
  img1.style.marginTop = "10px";
  card.appendChild(img1);

  // contact img
  const subtitle = document.createElement("div");
  subtitle.className = "muted heatmap-subtitle";
  subtitle.style.marginTop = "12px";
  subtitle.textContent = "Contact map (<= 8Å)";
  card.appendChild(subtitle);

  
}


function ssToClass(ss) {
  if (ss === "H") return "ss-H";
  if (ss === "E") return "ss-E";
  if (ss === "T") return "ss-T";
  if (ss === "S") return "ss-S";
  return "ss-dot";
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
    });
  });
}

// Pozovi setupTabs kad se DOM učita
document.addEventListener("DOMContentLoaded", setupTabs);

