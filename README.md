# Protein App

Django app for protein distance/analysis with a Neo4j backend and a small SQLite database.

## Requirements

- Python 3.12+
- Neo4j database with imported protein data

## Local quick start

1) Create a virtualenv and install deps:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2) Set env vars (example):

```bash
export NEO4J_URI=bolt://localhost:7687
export NEO4J_USER=neo4j
export NEO4J_PASS=your_password
# Optional
export DJANGO_SECRET_KEY=dev-only-unsafe-secret-key
export DJANGO_DEBUG=1
```

3) Run the app:

```bash
python manage.py migrate
python manage.py runserver
```

Open:

- App: `http://127.0.0.1:8000`

## Docker quick start

This project can run without Neo4j Desktop if you provide a dump of an existing Neo4j database.

### Expected structure

Put the dump file here:

```text
neo4j/backup/neo4j.dump
```

### 1. Prepare `.env`

Copy `.env.example` to `.env` and set values:

```env
NEO4J_USER=neo4j
NEO4J_PASS=your-password
DJANGO_SECRET_KEY=dev-only-unsafe-secret-key
DJANGO_DEBUG=1
DJANGO_ALLOWED_HOSTS=localhost,127.0.0.1
```

### 2. Start Neo4j only

```bash
docker compose up -d neo4j
```

### 3. Stop Neo4j before restore

```bash
docker compose stop neo4j
```

### 4. Restore the dump

If the database name is `neo4j`:

```bash
docker compose run --rm neo4j neo4j-admin database load neo4j --from-path=/backup --overwrite-destination=true
```

### 5. Start all services

```bash
docker compose up -d
```


Open:

- App: `http://localhost:8000`
- Neo4j Browser: `http://localhost:7474`

## Project notes

- SQLite is used for Django’s default database at `db.sqlite3`.
- Neo4j connection is required for query endpoints.
- The application reads Neo4j connection settings from `NEO4J_URI`, `NEO4J_USER`, and `NEO4J_PASS`.
- Docker setup is intended to use an existing Neo4j dump, not to re-import CSV files.
- For Docker restore details, see `DOCKER_NEO4J_DUMP.md`.
