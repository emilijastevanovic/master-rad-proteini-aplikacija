from django.urls import path
from . import views 

# putanja, hendler, ime
urlpatterns = [
    path("graph3d/", views.graph3d_protein, name="api_graph3d"),
    path("stats/", views.stats, name="api_stats"),
    path("global_stats/", views.global_stats, name="api_global_stats"),
    path("heatmap/distance/", views.heatmap_distance, name="api_heatmap_distance"),
    path("heatmap/segments/", views.segments_histogram, name="api_segments_histogram"),
]