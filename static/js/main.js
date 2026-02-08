let SEQ_CACHE = null;     
let SEQ_RENDERED = false;

let HEATMAP_URLS = null;

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
    if (!nodes.length) {
      alert("Nema rezultata za izabrane filtere (protein/lanac/prag).");
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
    if (e.key === "Escape") hideNodeTooltip();
  });
});


