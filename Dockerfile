FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# matplotlib (Agg backend) i neo4j driver trebaju ove pakete
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir gunicorn whitenoise

COPY . .

# collectstatic ne koristi Neo4j, ali Django treba SECRET_KEY
RUN DJANGO_SECRET_KEY=build-only-key \
    NEO4J_URI=bolt://localhost:7687 \
    NEO4J_USER=neo4j \
    NEO4J_PASS=build \
    python manage.py collectstatic --noinput

EXPOSE 8000

CMD ["gunicorn", "proteinapp.wsgi:application", \
     "--bind", "0.0.0.0:8000", \
     "--workers", "2", \
     "--timeout", "120"]
