from unittest.mock import Mock, patch

from django.core.cache import cache
from django.test import TestCase

from .query_builder import build_graph3d_protein
from .stats_queries import stat_aa_seq
from .views import STAT_CACHE_PARAMS, STAT_HANDLERS


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


class AaSeqTests(TestCase):
    ROWS = [
        {"name": "MET", "ss": "H"},
        {"name": "GLN", "ss": "H"},
        {"name": "ILE", "ss": "E"},
    ]

    @patch("distances.stats_queries.run_query")
    def test_sequences_come_from_single_aligned_query(self, run_query):
        run_query.return_value = self.ROWS

        result = stat_aa_seq("1a17", None, None, None, None)

        # jedan upit — front spaja aa_seq i aa_ss_seq po poziciji, pa oba
        # niza moraju poticati iz istog skupa redova
        self.assertEqual(run_query.call_count, 1)
        self.assertEqual(
            [r["name"] for r in result["aa_seq"]],
            ["MET", "GLN", "ILE"],
        )
        self.assertEqual(
            [r["ss"] for r in result["aa_ss_seq"]],
            ["H", "H", "E"],
        )
        self.assertEqual(len(result["aa_seq"]), len(result["aa_ss_seq"]))

    @patch("distances.stats_queries.run_query")
    def test_orders_by_chain_then_index(self, run_query):
        run_query.return_value = []

        stat_aa_seq("1a17", None, None, None, None)

        cypher = run_query.call_args.args[0]
        self.assertIn("ORDER BY aa.chain, aa.index", cypher)

    @patch("distances.stats_queries.run_query")
    def test_normalizes_protein_and_blank_chain(self, run_query):
        run_query.return_value = []

        result = stat_aa_seq(" 1a17 ", "", None, None, None)

        self.assertEqual(result["protein"], "1A17")
        self.assertEqual(run_query.call_args.kwargs["protein"], "1A17")
        self.assertIsNone(run_query.call_args.kwargs["chain"])


class StatsCacheTests(TestCase):
    def setUp(self):
        cache.clear()
        self.addCleanup(cache.clear)

    def _handler(self):
        return Mock(return_value={"total": 3, "aa_composition": [{"name": "ALA"}]})

    def _get(self, **params):
        query = {"protein": "1D3B", "include": "aa_composition"}
        query.update(params)
        return self.client.get("/api/stats/", query)

    def test_every_handler_declares_its_cache_key(self):
        # nova statistika mora svesno da odluči po čemu se kešira
        self.assertEqual(set(STAT_HANDLERS), set(STAT_CACHE_PARAMS))

    def test_repeated_request_is_served_from_cache(self):
        handler = self._handler()
        with patch.dict(STAT_HANDLERS, {"aa_composition": handler}):
            first = self._get(chain="A", type="caca")
            second = self._get(chain="A", type="caca")

        self.assertEqual(handler.call_count, 1)
        self.assertEqual(first.json()["aa_composition"], second.json()["aa_composition"])

    def test_different_chain_is_a_different_key(self):
        handler = self._handler()
        with patch.dict(STAT_HANDLERS, {"aa_composition": handler}):
            self._get(chain="A")
            self._get(chain="B")

        self.assertEqual(handler.call_count, 2)

    def test_distance_bounds_do_not_affect_key(self):
        # aa_composition ih ne koristi u upitu, pa ne smeju da dele keš
        handler = self._handler()
        with patch.dict(STAT_HANDLERS, {"aa_composition": handler}):
            self._get(chain="A", max_distance="8", min_distance="4")
            self._get(chain="A", max_distance="12", min_distance="")

        self.assertEqual(handler.call_count, 1)

    def test_type_affects_only_type_sensitive_stats(self):
        aa = self._handler()
        summary = Mock(return_value={"type": "caca", "count": 10})

        with patch.dict(STAT_HANDLERS, {"aa_composition": aa, "distance_summary": summary}):
            self._get(include="aa_composition,distance_summary", type="caca")
            self._get(include="aa_composition,distance_summary", type="minbezh")

        self.assertEqual(aa.call_count, 1)       # ne zavisi od tipa
        self.assertEqual(summary.call_count, 2)  # zavisi od tipa

    def test_blank_and_prazno_type_share_one_entry(self):
        summary = Mock(return_value={"type": "caca", "count": 10})
        with patch.dict(STAT_HANDLERS, {"distance_summary": summary}):
            self._get(include="distance_summary", type="")
            self._get(include="distance_summary", type="prazno")
            self._get(include="distance_summary", type="caca")

        self.assertEqual(summary.call_count, 1)

    def test_empty_result_is_not_cached(self):
        handler = Mock(return_value={"type": "caca", "empty": True})
        with patch.dict(STAT_HANDLERS, {"distance_summary": handler}):
            self._get(include="distance_summary")
            self._get(include="distance_summary")

        self.assertEqual(handler.call_count, 2)

    def test_unknown_stat_still_reported(self):
        response = self._get(include="nepostojeca")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["nepostojeca"], {"error": "unknown stat"})
