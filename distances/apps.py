from django.apps import AppConfig
import atexit

class DistancesConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "distances"

    def ready(self):
        from .neo4j_client import close_driver
        atexit.register(close_driver)
