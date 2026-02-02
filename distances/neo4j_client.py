import os
from neo4j import GraphDatabase
from dotenv import load_dotenv

load_dotenv()

_DRIVER = None

def _get_neo4j_config():
    uri = os.getenv("NEO4J_URI")
    user = os.getenv("NEO4J_USER")
    password = os.getenv("NEO4J_PASS")  # ili NEO4J_PASSWORD ako tako zoveš u .env
    missing = [k for k, v in {
        "NEO4J_URI": uri,
        "NEO4J_USER": user,
        "NEO4J_PASS": password,
    }.items() if not v]
    return uri, user, password, missing

def get_driver():
    global _DRIVER
    if _DRIVER is None:
        uri, user, password, missing = _get_neo4j_config()
        if missing:
            raise RuntimeError(f"Neo4j config missing env vars: {', '.join(missing)}")

        _DRIVER = GraphDatabase.driver(
            uri,
            auth=(user, password),
            max_connection_pool_size=50,
            connection_acquisition_timeout=10.0,
            connection_timeout=5.0,
            max_transaction_retry_time=10.0,
        )
    return _DRIVER

def close_driver():
    global _DRIVER
    if _DRIVER is not None:
        _DRIVER.close()
        _DRIVER = None

def run_query(query, **params):
    driver = get_driver()
    with driver.session(default_access_mode="READ") as session:
        result = session.run(query, parameters=params, timeout=10)
        return [r.data() for r in result]
