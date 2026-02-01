import os
from neo4j import GraphDatabase
from dotenv import load_dotenv

load_dotenv()

uri = os.getenv("NEO4J_URI")
user = os.getenv("NEO4J_USER")
password = os.getenv("NEO4J_PASS")

driver = GraphDatabase.driver(uri, auth=(user, password))

def run_query(query, **params):
    with driver.session() as session:
        result = session.run(query, **params)
        return [r.data() for r in result]