# Arhitektura projekta

Pregled strukture aplikacije — čemu služi svaki folder i fajl, i kako podaci
teku od klika u browseru do Neo4j baze i nazad.

Aplikacija je Django veb-aplikacija za 3D vizuelizaciju i statističku analizu
rastojanja između aminokiselina u proteinima. Podaci se čuvaju u Neo4j graf
bazi (čvorovi `AminoAcid`, veze `DISTANCE`), a Django služi kao sloj koji nad
njima izvršava Cypher upite i generiše slike.

## Tok jednog zahteva

```
browser  →  templates/index.html  +  static/js/main.js
              ↓ fetch /api/...
           proteinapp/urls.py  →  distances/urls.py  →  distances/views.py
                                                          ↓
                          distances/query_builder.py    (gradi Cypher za graf)
                          distances/stats_queries.py    (Cypher za statistike)
                                                          ↓
                          distances/neo4j_client.py     (drajver, izvršava upit)
                                                          ↓
                                                     Neo4j baza
                          distances/protein_analysis.py (pandas + matplotlib → PNG)
```

## Folderi

| Folder | Uloga |
|---|---|
| `proteinapp/` | Django *projekat* — globalna podešavanja i rutiranje |
| `distances/` | Django *aplikacija* — sva logika: pogledi, upiti, analiza |
| `templates/` | HTML šablon (samo jedan — aplikacija je jednostrana) |
| `static/` | CSS i JavaScript koje servira whitenoise |
| `neo4j/backup/` | `neo4j.dump` — bind-mount u kontejner na `/backup`, za `neo4j-admin restore` |
| `latex_heatmaps/`, `latex_segment_histograms/` | slike koje aplikacija **sama snima** dok radi, za 4 odabrana proteina |
| `overleaf_images/`, `overleaf_images_simple/` | ručno pripremljene slike za rad (+ `.zip` za upload na Overleaf) |
| `distances/migrations/` | prazan — Neo4j nije Django ORM, migracije postoje samo za ugrađene Django tabele |

## Backend — `proteinapp/` (projekat)

| Fajl | Šta radi |
|---|---|
| `settings.py` | 116 linija. Sve iz env promenljivih preko `python-dotenv`. Tu je `CACHES` (LocMemCache), `whitenoise` u MIDDLEWARE, `TEMPLATES.DIRS → templates/`, `INSTALLED_APPS` sa `distances.apps.DistancesConfig` |
| `urls.py` | Tri rute: `/` → `home()` (renderuje `index.html`), `/admin/`, `/api/` → uključuje `distances/urls.py` |
| `wsgi.py` | Ulazna tačka za gunicorn — ovo pokreće `docker-compose.yml` |
| `asgi.py` | Async varijanta, nekorišćena |
| `manage.py` | Standardni Django CLI (`migrate`, `test`, `collectstatic`) |

## Backend — `distances/` (aplikacija)

| Fajl | Linija | Šta radi |
|---|---|---|
| `views.py` | 422 | **Srce backenda.** 5 endpointa: `graph3d_protein` (podaci za 3D), `stats` (dispečer za 6 statistika preko `STAT_HANDLERS`), `global_stats` (4 upita nad celom bazom), `heatmap_distance` i `segments_histogram` (vraćaju PNG). Tu je i keširanje, plus `save_latex_*` pomoćne koje snimaju PNG na disk |
| `stats_queries.py` | 278 | **Svi Cypher upiti za statistike** — po jedna funkcija za svaku (`stat_aa_composition`, `stat_ss_distribution`, `stat_aa_seq`, `stat_ss_pairs`, `stat_distance_summary`, `stat_outlier_ss`, `stat_segment_lengths`, `stat_aa_heatmap_df`) |
| `query_builder.py` | 104 | **Samo za 3D graf.** Sklapa Cypher dinamički, jer graf ima 8 opcionih filtera. Dve grane: „node-only" kad nema nijednog filtera veza, i „edge grana" inače |
| `neo4j_client.py` | 47 | Drajver ka bazi. `get_driver()` pravi jedan globalni drajver lenjo (pool 50 konekcija), `run_query()` otvara READ sesiju i vraća listu rečnika |
| `protein_analysis.py` | 406 | **pandas + matplotlib.** Za aplikaciju su bitne samo 3 funkcije: `make_heat_maps` (matrica → PNG), `make_segment_hist_png` (histogram → PNG), `make_heatmap_matrix`. Ostalih ~14 su iz ranije, sveska-orijentisane faze rada |
| `apps.py` | 44 | **Pre-zagrevanje.** Pri startu pokreće pozadinsku nit koja unapred izračuna `global_stats`, sa 6 pokušaja na po 10s (čeka da Neo4j postane spreman). Takođe registruje zatvaranje drajvera na izlazu |
| `tests.py` | 156 | 13 testova: `Graph3DProteinTests`, `AaSeqTests`, `StatsCacheTests` |
| `urls.py` | 10 | Mapiranje 5 API putanja na poglede |
| `models.py`, `admin.py` | 1 | Prazni — nema Django modela, podaci su u Neo4j |

## Frontend

| Fajl | Linija | Šta radi |
|---|---|---|
| `templates/index.html` | 310 | Kontrole (filteri) na vrhu, pa `.layout` grid: levo 3D prikaz, desno `#statsPanel` sa karticom globalnih statistika i 7 tabova. Učitava `3Dmol-min.js` sa CDN-a i `main.js` |
| `static/js/main.js` | 1301 | 25 funkcija. Grubo: `loadGraph3D` (dohvat + orkestracija svih ostalih), `renderGraph` (crta sfere i cilindre u 3Dmol), `showNodeTab` (fokus analiza u browseru), po jedan `load*` za svaku statistiku, `renderSeqIfReady`/`renderHeatmapsFromUrls`/`renderSegmentsImg` (lenjo iscrtavanje po tabovima), `setupTabs`, `aaToClass`/`ssToClass` (bojenje) |
| `static/css/style.css` | 900 | 11 sekcija — kontrole, layout, tabovi, kartica globalnih statistika, fokus analiza, tooltip, lightbox, AA čipovi, responsive |

## API endpointi

| Putanja | Pogled | Vraća |
|---|---|---|
| `/api/graph3d/` | `graph3d_protein` | JSON — čvorovi i veze za 3D prikaz |
| `/api/stats/` | `stats` | JSON — jedna ili više statistika, bira se preko `include=` |
| `/api/global_stats/` | `global_stats` | JSON — statistike cele baze |
| `/api/heatmap/distance/` | `heatmap_distance` | PNG — toplotna mapa rastojanja |
| `/api/heatmap/segments/` | `segments_histogram` | PNG — histogram dužina segmenata |

## Konfiguracija i dokumentacija

| Fajl | Šta |
|---|---|
| `Dockerfile` | Python 3.12-slim, `COPY . .`, `collectstatic` pri build-u — otud potreba za rebuild-om posle izmene koda |
| `docker-compose.yml` | Dva servisa: `neo4j` (healthcheck) i `web` (čeka da neo4j bude zdrav, pa `migrate` + gunicorn) |
| `requirements.txt` | 33 paketa — Django 5.2.7, neo4j drajver, pandas, matplotlib, python-dotenv |
| `README.md` | Uputstvo za pokretanje i korišćenje |
| `DOCKER_NEO4J_DUMP.md` | Kako uvesti `neo4j.dump` bez Neo4j Desktop-a |
| `STATUS_APLIKACIJE.md` | Opis trenutnog stanja aplikacije |
| `IMPROVEMENTS.txt` | Stara lista uočenih problema |

## Podaci — šta preživljava restart

| Šta | Gde | Preživljava restart kontejnera |
|---|---|---|
| Neo4j baza | imenovani volume `neo4j_data` | ✓ |
| `neo4j.dump` | bind-mount `./neo4j/backup` | ✓ |
| LaTeX slike | bind-mount `./latex_*` | ✓ |
| Kod, `staticfiles/` | u Docker image-u | ✗ (pravi se pri build-u) |
| Keš (statistike, PNG) | RAM `gunicorn` procesa | ✗ |

## Napomene o stanju koda

**`IMPROVEMENTS.txt` je uglavnom zastareo.** Od 8 stavki, 6 je rešeno: `chain`
filter u `query_builder.py` sada radi, `SECRET_KEY`/`DEBUG` idu iz env-a,
drajver se pravi lenjo, `graph3d` validira protein, `make_aa_string` je izvučen
na nivo modula, `loadGraph()` više ne postoji. Ostaje da su
`django-cors-headers` i `djangorestframework` u `requirements.txt` a nisu u
`INSTALLED_APPS` — dakle instaliraju se bez potrebe.

**`protein_analysis.py` je najveći fajl a najmanje se koristi** — 3 od 17
funkcija su povezane sa aplikacijom. Ostale (`make_distance_statistic_df`,
`topn_with_other`, `print_ss_tracks`, `plot_pie_from_percent_table`…) su iz
faze kad se analiza radila u svesci, pre nego što je logika prešla u Cypher.

**Statistike ne poštuju filtere rastojanja.** `max_distance` i `min_distance`
utiču samo na 3D prikaz. Statistike u bočnom panelu se računaju nad celim
proteinom (uz poštovanje lanca i tipa rastojanja). Isto važi za filtere `ss`,
`aminoname` i `k`, koji deluju samo na graf.

**Tab „Sekvenca" namerno prikazuje ceo protein**, bez obzira na izabrani lanac.
