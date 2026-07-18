from unittest.mock import patch

from django.test import TestCase

from .query_builder import build_graph3d_protein


class Graph3DProteinTests(TestCase):
    @patch("distances.views.run_query", return_value=[])
    def test_view_normalizes_protein_id(self, run_query):
        response = self.client.get("/api/graph3d/", {"protein": " 1a17 "})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(run_query.call_args.kwargs["protein"], "1A17")

    def test_query_filters_requested_protein(self):
        query, params = build_graph3d_protein({
            "protein": "2XYZ",
            "aminonames": [],
            "ss": "prazno",
            "k": "",
            "type": "prazno",
            "max_distance": "",
            "min_distance": "",
            "chain": "prazno",
        })

        self.assertIn("a.protein = $protein", query)
        self.assertEqual(params["protein"], "2XYZ")
