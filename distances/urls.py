from django.urls import path
from . import views 

urlpatterns = [
    path("graph3d/", views.graph3d_protein, name="api_graph3d"),
    path("chains/", views.chains, name="api_chains"),
    path("stats/", views.stats, name="api_stats"),
    path("global_stats/", views.global_stats, name="api_global_stats"),
    path("heatmap/distance/", views.heatmap_distance, name="api_heatmap_distance"),
    path("heatmap/aa_types/", views.heatmap_aa_types, name="api_heatmap_aa_types"),
    path("heatmap/segments/", views.segments_histogram, name="api_segments_histogram"),
    path("heatmap/segments_by_ss/", views.segments_histogram_by_ss, name="api_segments_histogram_by_ss"),
]
