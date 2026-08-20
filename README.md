# Protein App

Veb-aplikacija za 3D vizuelizaciju i statističku analizu aminokiselina i
njihovih rastojanja u proteinima. Aplikacija je napravljena u Django-u, dok se
podaci o proteinima čuvaju u Neo4j bazi.

## Preduslovi

Preporučeni način pokretanja je pomoću Docker-a. Potrebni su:

- Docker i Docker Compose;
- Neo4j dump baze `neo4j.dump`, koji se zbog veličine dostavlja odvojeno od
  izvornog koda.

Za lokalno pokretanje bez Docker-a potrebni su Python 3.12+ i već pokrenuta
Neo4j instanca u koju su uvezeni podaci.

## Baza podataka

Dump baze nije deo ovog repozitorijuma jer zauzima oko 2 GB. Objavljen je u
zasebnom repozitorijumu
[protein-distance-graph-db](https://github.com/emilijastevanovic/protein-distance-graph-db)
i preuzima se ovako:

```bash
curl -L -o neo4j/backup/neo4j.dump \
  https://github.com/emilijastevanovic/protein-distance-graph-db/releases/latest/download/neo4j.dump
```

Dump treba da se nađe na sledećoj putanji u projektu:

```text
neo4j/backup/neo4j.dump
```

Bez ovog dump-a moguće je pokrenuti Django server, ali prikaz proteina i
statistički upiti neće raditi.

## Pokretanje pomoću Docker-a

1. Kopirati primer konfiguracije:

```bash
cp .env.example .env
```

Po potrebi promeniti lozinku i ostale vrednosti u `.env` fajlu. Vrednost
`NEO4J_PASS` mora biti ista za Neo4j i Django servis.

2. Proveriti da se dump nalazi na putanji:

```text
neo4j/backup/neo4j.dump
```

3. Kreirati Neo4j volume i zatim zaustaviti servis pre učitavanja dump-a:

```bash
docker compose up -d neo4j
docker compose stop neo4j
```

4. Učitati dump u bazu:

```bash
docker compose run --rm neo4j neo4j-admin database load neo4j --from-path=/backup --overwrite-destination=true
```

Ovaj korak je potreban samo prilikom prvog pokretanja sa praznim Docker
volume-om.

5. Pokrenuti aplikaciju i bazu:

```bash
docker compose up -d
```

6. Otvoriti:

- aplikaciju: `http://localhost:8000`;
- Neo4j Browser: `http://localhost:7474`.

Status servisa i njihove logove moguće je proveriti komandama:

```bash
docker compose ps
docker compose logs -f web
docker compose logs -f neo4j
```

Aplikacija se zaustavlja komandom:

```bash
docker compose down
```

Komanda `docker compose down -v` briše i Neo4j volume i koristi se samo kada je
potrebno ponovo učitati dump od početka.

## Lokalno pokretanje bez Docker-a

1. Kreirati i aktivirati virtuelno okruženje, pa instalirati zavisnosti:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Na Windows-u se virtuelno okruženje aktivira komandom:

```powershell
.venv\Scripts\activate
```

2. Kopirati konfiguraciju i prilagoditi adresu i pristupne podatke postojećoj
Neo4j instanci:

```bash
cp .env.example .env
```

Za lokalnu Neo4j instancu podrazumevana adresa je `bolt://localhost:7687`.

3. Pokrenuti migracije i razvojni server:

```bash
python manage.py migrate
python manage.py runserver
```

Aplikacija je dostupna na `http://127.0.0.1:8000`.

## Korišćenje aplikacije

1. U polje **ID proteina** uneti identifikator proteina, na primer `1A17`.
2. Po želji izabrati aminokiseline i postaviti filtere:
   tip rastojanja, minimalno i maksimalno rastojanje, sekundarnu strukturu,
   minimalnu sekvencijalnu udaljenost `k` i lanac.
3. Kliknuti na **Prikaži**. Levo se prikazuje 3D graf, a desno statistike.
4. Prelaskom preko čvora prikazuju se osnovni podaci o aminokiselini. Klikom na
   čvor otvara se tab **Čvor**, sa kontaktima i analizom njegovog susedstva.
5. Ostali tabovi prikazuju:
   sastav aminokiselina i sekundarnih struktura, sekvence, pregled i toplotnu
   mapu rastojanja, statistike parova sekundarnih struktura, dužine segmenata i
   odstupajuće rezidue.
6. Klikom na toplotnu mapu ili histogram slika se otvara u uvećanom prikazu.

Ako kombinacija filtera nema rezultate, aplikacija prikazuje odgovarajuću
poruku. Dugme **Sve** bira sve aminokiseline, a **Ništa** poništava njihov izbor.

## Primeri za proveru

Dump korišćen tokom razvoja sadrži, između ostalih, sledeće proteinske ID-jeve:

- `1A17`
- `1K3I`
- `1G2Y`
- `1FEW`

Najjednostavnija provera je unos `1A17` bez dodatnih filtera i klik na
**Prikaži**. Za proveru relacija može se izabrati tip rastojanja `Cα–Cα` i
maksimalno rastojanje `8.0 Å`.

## Provera projekta

Posle instalacije zavisnosti mogu se izvršiti Django provere i testovi:

```bash
python manage.py check
python manage.py test
```

SQLite se koristi za podrazumevane Django tabele u lokalnom fajlu `db.sqlite3`.
Neo4j konekcija je neophodna za podatke o proteinima, graf i statističke API
rute. Dodatni detalji o učitavanju dump-a nalaze se u
[`DOCKER_NEO4J_DUMP.md`](DOCKER_NEO4J_DUMP.md).
