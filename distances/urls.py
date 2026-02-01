from django.urls import path
from . import views 

# putanja, hendler, ime
urlpatterns = [
    # API
    path("graph3d/", views.graph3d_protein, name="api_graph3d"),
    path("stats/", views.stats, name="api_stats"),
    path("heatmap/distance/", views.heatmap_distance, name="api_heatmap_distance"),
]