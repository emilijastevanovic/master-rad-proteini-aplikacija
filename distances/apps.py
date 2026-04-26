from django.apps import AppConfig
import atexit
import threading
import logging

logger = logging.getLogger(__name__)

class DistancesConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "distances"

    def ready(self):
        from .neo4j_client import close_driver
        atexit.register(close_driver)
        threading.Thread(target=_prewarm_global_stats, daemon=True).start()


def _prewarm_global_stats():
    import time
    # čeka da Neo4j bude spreman (do 60s)
    for attempt in range(6):
        try:
            from django.core.cache import cache
            if cache.get("global_stats"):
                return
            from distances.views import global_stats
            from django.test import RequestFactory
            global_stats(RequestFactory().get("/api/global_stats/"))
            logger.info("global_stats pre-warm done")
            return
        except Exception as e:
            logger.warning("global_stats pre-warm attempt %d failed: %s", attempt + 1, e)
            time.sleep(10)
