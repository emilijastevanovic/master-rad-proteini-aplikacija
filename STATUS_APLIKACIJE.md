# Status aplikacije

## Kratak opis

Ovo je Django aplikacija za pregled i analizu proteinskih podataka nad Neo4j graf bazom.  
Frontend prikazuje 3D graf aminokiselina i njihove distance, a backend vraća statistike i generiše slike (heatmap i histogram).

Trenutno aplikacija koristi:

- `Django` za web aplikaciju
- `Neo4j` kao glavnu bazu za proteinske podatke
- `SQLite` samo za Django internu bazu (`db.sqlite3`)
- `3Dmol.js` za 3D prikaz proteina
- `pandas`, `matplotlib`, `seaborn`, `numpy` za analitiku i vizuelizaciju

## Šta je implementirano

### 1. Osnovna web aplikacija

- Projekat je organizovan kao Django projekat `proteinapp`
- Glavna aplikacija je `distances`
- Root ruta `/` otvara stranicu `templates/index.html`
- API rute su pod `/api/`

## 2. Neo4j konekcija

Implementirana je konekcija ka Neo4j preko env varijabli:

- `NEO4J_URI`
- `NEO4J_USER`
- `NEO4J_PASS`

Konekcija je centralizovana u fajlu `distances/neo4j_client.py`.

Podržano je:

- validiranje da li postoje env varijable
- kreiranje Neo4j driver-a
- izvršavanje read-only Cypher upita
- zatvaranje driver-a pri gašenju aplikacije

## 3. 3D prikaz proteinskog grafa

Implementirano u:

- backend: `distances/views.py`, `distances/query_builder.py`
- frontend: `static/js/main.js`, `templates/index.html`

Funkcionalnosti:

- unos protein ID-ja
- filtriranje po aminokiselini
- filtriranje po tipu distance:
  - `caca`
  - `minbezh`
  - `maxbezh`
- filtriranje po lancu
- filtriranje po minimalnoj i maksimalnoj distanci
- prikaz čvorova aminokiselina u 3D prostoru
- prikaz veza između aminokiselina
- tooltip za čvorove:
  - protein
  - lanac
  - indeks
  - aminokiselina
  - sekundarna struktura
- tooltip za grane:
  - tip distance
  - vrednost distance

Napomena:

- ako nisu zadati filteri za distance, backend može da vrati samo čvorove sa koordinatama
- ako jesu zadati filteri, vraća i čvorove i veze

## 4. Statistike nad proteinom

API ruta: `/api/stats/`

Podržane statistike:

### AA composition

- broj pojavljivanja svake aminokiseline
- procenat svake aminokiseline u proteinu

### Sekvenca

- vraćanje primarne sekvence po redosledu indeksa
- vraćanje SS sekvence po redosledu indeksa

### SS pair statistika

Za parove sekundarnih struktura iz distance relacija računa se:

- `count`
- `mean`
- `median`
- `min`
- `max`
- `std`

## 5. Heatmap i histogram

Implementirano u `distances/views.py` i `distances/protein_analysis.py`.

Podržano:

- generisanje PNG heatmap slike matrice distanci
- generisanje histograma dužina segmenata sekundarne strukture
- cache rezultata preko Django cache-a

API rute:

- `/api/heatmap/distance/`
- `/api/heatmap/segments/`

Važno:

- kod za contact map postoji u analitičkom modulu, ali trenutno nije vraćen iz API-ja
- trenutno se ka frontend-u vraća samo `distance_png`

## 6. Frontend interfejs

Na stranici trenutno postoje tabovi:

- `Overview`
- `Seq`
- `Distance`
- `SS pairs`
- `Segments`

Frontend trenutno radi sledeće:

- šalje pozive ka backend API-jima
- učitava 3D graf
- učitava AA composition
- učitava SS pair statistiku
- priprema URL za heatmap
- priprema URL za histogram segmenata
- prefetch sekvence za prikaz u `Seq` tabu

## 7. Konfiguracija i deploy detalji

Urađeno je sledeće:

- `SECRET_KEY`, `DEBUG` i `ALLOWED_HOSTS` su prebačeni na env varijable
- statički fajlovi su podešeni preko `WhiteNoise`
- postoji README sa osnovnim koracima za pokretanje

## Šta model podataka verovatno sadrži u Neo4j

Na osnovu postojećih Cypher upita, aplikacija očekuje:

### Čvorove `AminoAcid`

Sa atributima približno ovog tipa:

- `protein`
- `chain`
- `index`
- `name`
- `ss`
- `ca_coordinates.x`
- `ca_coordinates.y`
- `ca_coordinates.z`

### Relacije `DISTANCE`

Sa atributima:

- `type`
- `value`

To znači da je Neo4j baza već morala biti pripremljena ranije, ali skripta za taj import trenutno nije vidljiva u ovom repozitorijumu.

## Gde si verovatno stala

Po stanju koda izgleda da je završen:

- osnovni Django projekat
- konektor prema Neo4j
- glavni 3D viewer
- osnovne statistike
- heatmap distance
- histogram segmenata
- chain filter kroz API i graf upite
- osnovno keširanje i konfiguracija

Deluje da nisu završeni ili nisu potpuno dovedeni do kraja:

- kompletan prikaz `Distance summary` u UI-ju (`TODO` još stoji u template-u)
- vraćanje i prikaz `contact map`
- testovi (`distances/tests.py` je prazan)
- Django modeli nisu realno korišćeni (`distances/models.py` je prazan)
- skripta/upit za inicijalni punjenje Neo4j baze nije pronađen
- CSV ulazni fajlovi nisu u ovom repozitorijumu

## Neo4j / CSV deo koji je bitan za nastavak rada

Važna napomena:

- u trenutnom repozitorijumu nisam našla nijedan `.csv` fajl
- nisam našla import skriptu za Neo4j
- nisam našla `LOAD CSV`, `.cypher`, `.cql` ili sličan import fajl
- u git istoriji ovog repozitorijuma se takođe ne vidi import skripta

To znači da je jedan od sledećih scenarija verovatan:

1. CSV fajlovi sa 933 proteina bili su van repozitorijuma
2. import u Neo4j je rađen ručno u Neo4j Browser-u
3. import skripta je bila u drugom folderu/projektu
4. baza je napunjena lokalno, ali kod za import nije sačuvan ovde

## Šta bi sledeće trebalo pronaći

Da bi mogla da nastaviš rad bez nagađanja, treba pronaći jedno od sledećeg:

- originalne CSV fajlove
- Cypher `LOAD CSV` upit
- Python skriptu koja je unosila podatke u Neo4j
- export Neo4j baze
- beleške ili stare fajlove sa mapiranjem CSV kolona na `AminoAcid` i `DISTANCE`

## Korisni fajlovi za nastavak

- `README.md`
- `distances/views.py`
- `distances/query_builder.py`
- `distances/stats_queries.py`
- `distances/protein_analysis.py`
- `distances/neo4j_client.py`
- `static/js/main.js`
- `templates/index.html`

## Kratak zaključak

Aplikacija nije početna skica, već funkcionalan prototip sa Neo4j upitima, 3D prikazom i osnovnim statističkim modulima.  
Najveća rupa za nastavak rada trenutno nije frontend ili Django sloj, nego nedostatak izvora za inicijalno punjenje Neo4j baze iz CSV podataka za 933 proteina.
