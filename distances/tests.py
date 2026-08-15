from unittest.mock import Mock, patch

from django.core.cache import cache
from django.test import TestCase

import pandas as pd

from .protein_analysis import make_segment_hist_png, make_aa_type_matrix, make_aa_type_heatmap_png
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


class DistanceParamValidationTests(TestCase):
    BASE = {
        "protein": "1G2Y", "aminoname": "", "chain": "",
        "ss": "prazno", "k": "", "type": "caca",
        "max_distance": "", "min_distance": "",
    }

    def _get(self, **params):
        query = dict(self.BASE)
        query.update(params)
        return self.client.get("/api/graph3d/", query)

    @patch("distances.views.run_query", return_value=[])
    def test_comma_is_accepted_as_decimal_separator(self, run_query):
        # na našoj tastaturi je „8,5" prirodnije od „8.5"
        response = self._get(max_distance="8,5")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(run_query.call_args.kwargs["max_distance"], 8.5)

    @patch("distances.views.run_query", return_value=[])
    def test_comma_accepted_for_lower_bound(self, run_query):
        response = self._get(min_distance="4,25")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(run_query.call_args.kwargs["min_distance"], 4.25)

    @patch("distances.views.run_query", return_value=[])
    def test_invalid_distance_returns_400_not_500(self, run_query):
        response = self._get(max_distance="abc")

        self.assertEqual(response.status_code, 400)
        self.assertIn("gornja granica", response.json()["detail"])
        run_query.assert_not_called()

    @patch("distances.views.run_query", return_value=[])
    def test_invalid_k_returns_400_not_500(self, run_query):
        response = self._get(k="8.5")

        self.assertEqual(response.status_code, 400)
        self.assertIn("sekvencijalno rastojanje", response.json()["detail"])
        run_query.assert_not_called()

    @patch("distances.views.run_query", return_value=[])
    def test_blank_bounds_are_not_applied(self, run_query):
        response = self._get(max_distance="", min_distance="prazno", k="")

        self.assertEqual(response.status_code, 200)
        self.assertNotIn("max_distance", run_query.call_args.kwargs)
        self.assertNotIn("min_distance", run_query.call_args.kwargs)
        self.assertNotIn("k", run_query.call_args.kwargs)


class SegmentHistogramTests(TestCase):
    def test_integer_bins_when_range_is_narrow(self):
        # opseg 1-8 ne sme da se deli na 30 pregrada
        png = make_segment_hist_png([1, 1, 2, 3, 5, 8], protein="TEST")

        self.assertTrue(png.startswith(b"\x89PNG"))

    def test_wide_range_still_renders(self):
        png = make_segment_hist_png(list(range(1, 200)), protein="TEST")

        self.assertTrue(png.startswith(b"\x89PNG"))


class AaTypeMatrixTests(TestCase):
    """Mapa proseka po vrsti aminokiseline (GLU × LYS × THR)."""

    def _df(self, rows):
        return pd.DataFrame(rows, columns=["aa1", "idx1", "aa2", "idx2", "distance"])

    def test_cell_is_mean_over_all_pairs_of_two_types(self):
        df = self._df([
            ("GLU", 1, "LYS", 2, 10.0),
            ("GLU", 3, "LYS", 4, 20.0),
        ])

        M = make_aa_type_matrix(df)

        self.assertAlmostEqual(M.loc["GLU", "LYS"], 15.0)

    def test_matrix_is_symmetric_although_query_returns_one_direction(self):
        # upit vraća par samo u jednom smeru; obe ćelije moraju biti popunjene
        df = self._df([("GLU", 1, "LYS", 2, 10.0)])

        M = make_aa_type_matrix(df)

        self.assertAlmostEqual(M.loc["GLU", "LYS"], M.loc["LYS", "GLU"])

    def test_diagonal_averages_pairs_of_same_type(self):
        df = self._df([
            ("THR", 1, "THR", 2, 4.0),
            ("THR", 2, "THR", 3, 6.0),
        ])

        M = make_aa_type_matrix(df)

        self.assertAlmostEqual(M.loc["THR", "THR"], 5.0)

    def test_axes_carry_the_same_names_in_the_same_order(self):
        df = self._df([
            ("THR", 1, "GLU", 2, 8.0),
            ("LYS", 3, "GLU", 4, 9.0),
        ])

        M = make_aa_type_matrix(df)

        self.assertEqual(list(M.index), list(M.columns))
        self.assertEqual(list(M.index), ["GLU", "LYS", "THR"])

    def test_png_is_rendered(self):
        df = self._df([("GLU", 1, "LYS", 2, 10.0)])

        png = make_aa_type_heatmap_png(df, protein="TEST")

        self.assertTrue(png.startswith(b"\x89PNG"))

    def test_empty_input_yields_no_image(self):
        png = make_aa_type_heatmap_png(self._df([]), protein="TEST")

        self.assertIsNone(png)


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
