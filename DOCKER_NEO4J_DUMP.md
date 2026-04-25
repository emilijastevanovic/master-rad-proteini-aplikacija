# Docker + Neo4j dump

Ovaj projekat moze da se pokrene bez Neo4j Desktop-a.

## Struktura

U projekat stavi dump baze ovde:

```text
neo4j/backup/neo4j.dump
```

Ako je ime baze drugacije, naziv dump fajla ce pratiti to ime.

## 1. Priprema `.env`

Kopiraj `.env.example` u `.env` i podesi vrednosti:

```env
NEO4J_USER=neo4j
NEO4J_PASS=your-password
DJANGO_SECRET_KEY=dev-only-unsafe-secret-key
DJANGO_DEBUG=1
DJANGO_ALLOWED_HOSTS=localhost,127.0.0.1
```

## 2. Podigni samo Neo4j

```bash
docker compose up -d neo4j
```

## 3. Zaustavi Neo4j pre restore-a

```bash
docker compose stop neo4j
```

## 4. Restore dump-a

Ako je ime baze `neo4j` i dump fajl je `neo4j.dump`:

```bash
docker compose run --rm neo4j neo4j-admin database load neo4j --from-path=/backup --overwrite-destination=true
```

Ako koristis drugo ime baze, zameni `neo4j` tim imenom.

## 5. Podigni sve servise

```bash
docker compose up -d
```

## 6. Provera

- Django app: `http://localhost:8000`
- Neo4j Browser: `http://localhost:7474`

## Napomene

- Ovo ne radi ponovni import iz CSV fajlova.
- Koristi se postojeci dump tvoje baze.
- Neo4j Desktop nije potreban na drugom racunaru.
- Prvi restore se radi samo jednom po praznom volume-u.
- Ako zelis ponovni restore na istom racunaru, prvo obrisi volume:

```bash
docker compose down -v
```
