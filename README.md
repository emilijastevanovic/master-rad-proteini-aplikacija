# Protein App

Django app for protein distance/analysis with a Neo4j backend and a small SQLite database.

## Quick start

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

## Notes

- SQLite is used for Django’s default database at `db.sqlite3`.
- Neo4j connection is required for query endpoints.
