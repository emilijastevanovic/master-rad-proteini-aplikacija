import logging

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io

logger = logging.getLogger(__name__)


def make_residual_level_df(df:pd.DataFrame) -> pd.DataFrame:

  """  
  Funkcija vraca DataFrame gde svaki red predstavlja pojedinacnu aminokiselinu u sekvenci proteina.
  """
  df = df.copy()

  res1 = (
    df[["protein", "lanacaa1", "aa1", "indeksaa1", "ss1", "idss1", "sclass1"]]
    .rename(columns={
        "lanacaa1": "lanac",
        "aa1": "aa",
        "indeksaa1": "indeks",
        "ss1": "ss",
        "idss1": "idss",
        "sclass1": "sclass",
    })
  )

  res2 = (
    df[["protein", "lanacaa2", "aa2", "indeksaa2", "ss2", "idss2", "sclass2"]]
    .rename(columns={
        "lanacaa2": "lanac",
        "aa2": "aa",
        "indeksaa2": "indeks",
        "ss2": "ss",
        "idss2": "idss",
        "sclass2": "sclass",
    })
  )

  res = pd.concat([res1, res2], ignore_index=True).drop_duplicates(
    subset=["protein", "lanac", "indeks"]
  ).reset_index(drop=True)

  res["node"] = res["lanac"].astype(str) + ":" + res["aa"].astype(str) + res["indeks"].astype(str)


  return res


def make_heatmap_matrix(df:pd.DataFrame) -> pd.DataFrame:
  """
  Funkcija pravi heatmap matricu rastojanja izmedju aminokiselina.
  """
  df_pom = df[["aa1", "idx1", "aa2", "idx2", "distance"]]
  df_pom["node1"] = df_pom["aa1"].astype(str) + df_pom["idx1"].astype(str)
  df_pom["node2"] = df_pom["aa2"].astype(str) + df_pom["idx2"].astype(str)

  matrix = df_pom.pivot_table(index="node1", columns="node2", values="distance", aggfunc="mean")
  matrix = matrix.combine_first(matrix.T)
  np.fill_diagonal(matrix.values, 0.0)

  return matrix

def make_distance_statistic_df(
    df: pd.DataFrame,
    distance_col: str = "rastojanje",
    percentiles=(1, 5, 10, 25, 50, 75, 90, 95, 99),
    thresholds=(4.0, 5.0, 8.0),
    add_outliers_iqr: bool = True
) -> pd.Series:

    x = pd.to_numeric(df[distance_col], errors="coerce").dropna().to_numpy()
    n_aa = len(make_residual_level_df(df))
    n = len(x)

    def _mad(arr: np.ndarray) -> float:
        med = np.median(arr)
        return np.median(np.abs(arr - med))

    if n == 0:
        out = {"count": 0}
        for k in ["mean","std","min","q25","median","q75","iqr","max","mad"]:
            out[k] = np.nan
        for p in percentiles:
            out[f"p{int(p):02d}"] = np.nan
        for t in thresholds:
            out[f"pct_lt_{t:g}"] = np.nan
        if add_outliers_iqr:
            out["n_outliers_iqr"] = np.nan
            out["pct_outliers_iqr"] = np.nan
        return pd.Series(out)

    q25 = float(np.percentile(x, 25))
    q50 = float(np.percentile(x, 50))
    q75 = float(np.percentile(x, 75))
    iqr = q75 - q25

    out = {
        "count": n,
        "count aa": n_aa,
        "mean": float(np.mean(x)),
        "std": float(np.std(x, ddof=1)) if n > 1 else 0.0,
        "min": float(np.min(x)),
        "q25": q25,
        "median": q50,
        "q75": q75,
        "iqr": float(iqr),
        "max": float(np.max(x)),
        "mad": float(_mad(x)),
    }

    for p in percentiles:
        out[f"p{int(p):02d}"] = float(np.percentile(x, p))

    for t in thresholds:
        out[f"pct_lt_{t:g}"] = float(np.mean(x < t)) * 100.0

    if add_outliers_iqr:
        lower = q25 - 1.5 * iqr
        upper = q75 + 1.5 * iqr
        n_out = int(np.sum((x < lower) | (x > upper)))
        out["n_outliers_iqr"] = n_out
        out["pct_outliers_iqr"] = (n_out / n) * 100.0

    return pd.Series(out)


def make_percent_table(df: pd.DataFrame, col: str) -> pd.DataFrame:
    t = (
        df[col]
        .fillna("NA")
        .value_counts(dropna=False)
        .rename_axis(col)
        .reset_index(name="count")
    )
    t["pct"] = 100 * t["count"] / t["count"].sum()
    return t

def topn_with_other(df: pd.DataFrame, col: str, n: int = 10) -> pd.DataFrame:
    t = make_percent_table(df, col)
    if len(t) <= n:
        return t
    top = t.iloc[:n].copy()
    other_count = t.iloc[n:]["count"].sum()
    other_pct = t.iloc[n:]["pct"].sum()
    top = pd.concat([top, pd.DataFrame([{col: "Other", "count": other_count, "pct": other_pct}])], ignore_index=True)
    return top


def count_segments_by_ss(res: pd.DataFrame) -> pd.DataFrame:
    seg = (
        res.dropna(subset=["idss"])
        .groupby(["protein", "lanac", "idss"])["ss"]
        .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else x.iloc[0])
        .reset_index(name="segment_ss")
    )

    out = (
        seg.groupby(["protein", "segment_ss"])
        .size()
        .reset_index(name="n_segments")
        .sort_values(["protein", "n_segments"], ascending=[True, False])
    )
    return out

def segment_length_table(res: pd.DataFrame) -> pd.DataFrame:
    seg_len = (
        res.dropna(subset=["idss"])
        .groupby(["protein", "lanac", "idss"])
        .size()
        .reset_index(name="segment_len")
        .sort_values(["protein", "segment_len"], ascending=[True, False])
    )
    return seg_len

def chain_summary(res: pd.DataFrame) -> pd.DataFrame:
    out = (
        res.groupby(["protein", "lanac"])
        .size()
        .reset_index(name="n_residues")
        .sort_values(["protein", "n_residues"], ascending=[True, False])
    )
    return out

def build_ss_string_with_gaps(chain_df: pd.DataFrame, ss_col="ss", index_col="indeks",
                              gap_char=".", fill_char="NA") -> str:
    """
    chain_df: df samo za jedan protein+lanac, sa kolonama ss i indeks
    Vraca string sekundarne strukture, uz ubacivanje gap_char gde fali indeks.
    """
    chain_df = chain_df.copy()
    chain_df[ss_col] = chain_df[ss_col].fillna(fill_char).astype(str)
    chain_df[index_col] = pd.to_numeric(chain_df[index_col], errors="coerce")

    chain_df = chain_df.dropna(subset=[index_col]).sort_values(index_col)
    if chain_df.empty:
        return ""

    idx = chain_df[index_col].astype(int).to_numpy()
    ss_vals = chain_df[ss_col].to_numpy()

    idx2ss = dict(zip(idx, ss_vals))

    start = idx.min()
    end = idx.max()

    out = []
    for i in range(start, end + 1):
        out.append(idx2ss.get(i, gap_char))

    return "".join(out)

def format_ss_track(ss_string: str,
                    start_index: int,
                    line_width: int = 80,
                    tick_every: int = 10,
                    label: str = "ss") -> str:
    """
    Pravi blok prikaz:
    0         10        20 ...
    ss: ---EEEHHH...
    """
    if not ss_string:
        return f"{label}: <empty>"

    blocks = []
    n = len(ss_string)

    for offset in range(0, n, line_width):
        chunk = ss_string[offset: offset + line_width]

        idx_line = [" "] * len(chunk)
        for pos in range(0, len(chunk), tick_every):
            num = str(start_index + offset + pos)
            for k, ch in enumerate(num):
                if pos + k < len(idx_line):
                    idx_line[pos + k] = ch
        idx_line = "".join(idx_line)

        ss_line = f"{label}: " + chunk

        blocks.append(idx_line)
        blocks.append(ss_line)
        blocks.append("")

    return "\n".join(blocks).rstrip()

def print_ss_tracks(res: pd.DataFrame,
                    ss_col="ss",
                    index_col="indeks",
                    protein_col="protein",
                    chain_col="lanac",
                    line_width: int = 80,
                    tick_every: int = 10,
                    gap_char: str = ".",
                    fill_char: str = "-",
                    max_chars: int | None = None):
    """
    Ispisuje SS trackove za svaki protein i lanac.
    gap_char: ubacuje se kad fali indeks (npr '.')
    fill_char: ako ss NaN, sta ubaciti (npr '-')
    max_chars: ako zelis skratiti prikaz (npr 300)
    """

    res = res.copy()

    res[index_col] = pd.to_numeric(res[index_col], errors="coerce")

    groups = res.dropna(subset=[protein_col, chain_col]).groupby([protein_col, chain_col], sort=True)

    for (prot, chain), g in groups:
        g = g.dropna(subset=[index_col]).sort_values(index_col)
        if g.empty:
            logger.debug("%s chain %s: empty", prot, chain)
            continue

        start = int(g[index_col].min())
        ss_string = build_ss_string_with_gaps(
            g,
            ss_col=ss_col,
            index_col=index_col,
            gap_char=gap_char,
            fill_char=fill_char
        )

        if max_chars is not None:
            ss_string = ss_string[:max_chars]

        track = format_ss_track(
            ss_string,
            start_index=start,
            line_width=line_width,
            tick_every=tick_every,
            label="ss"
        )

        logger.debug("%s chain %s  start_index=%d  length=%d\n%s", prot, chain, start, len(ss_string), track)

def plot_segment_length_hist(res: pd.DataFrame, bins: int = 30):
    seg_len = segment_length_table(res)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(seg_len["segment_len"], bins=bins)
    ax.set_title("Segment length distribution (residues per idss)")
    ax.set_xlabel("Duzina segmenta")
    ax.set_ylabel("Broj pojavljivanja")
    plt.tight_layout()
    plt.show()

def plot_pie_from_percent_table(t: pd.DataFrame, name_col: str, value_col: str = "count",
                                title: str = "", donut: bool = False):
    labels = t[name_col].astype(str).tolist()
    values = t[value_col].values

    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(
        values,
        labels=labels,
        autopct=lambda p: f"{p:.1f}%" if p > 2 else "",
        startangle=90
    )
    ax.set_title(title)

    if donut:
        centre_circle = plt.Circle((0, 0), 0.60, fc="white")
        ax.add_artist(centre_circle)

    ax.axis("equal")
    plt.tight_layout()
    plt.show()


def make_heat_maps(df: pd.DataFrame, protein: str):
  M = make_heatmap_matrix(df)

  import re
  def extract_num(s):
    m = re.search(r"(\d+)$", s)
    return int(m.group(1)) if m else 10**9

  M = M.loc[sorted(M.index, key=extract_num), sorted(M.columns, key=extract_num)]

  fig1, ax1 = plt.subplots(figsize=(12, 10))
  sns.heatmap(M, cmap="viridis", square=True, ax=ax1)
  ax1.set_title(f"Matrica rastojanja — {protein}")
  ax1.set_xlabel("Aminokiselina")
  ax1.set_ylabel("Aminokiselina")
  fig1.tight_layout()

  buf1 = io.BytesIO()
  fig1.savefig(buf1, format="png", dpi=150, bbox_inches="tight")
  plt.close(fig1)
  buf1.seek(0)


  return {
        "distance_png": buf1.getvalue(),
    }

def make_aa_string(df:pd.DataFrame) -> str:
    """
    Funkcija prima rezudial level DataFrame i vraca sekvencu aminokiselina.
    """
    return "-".join(list(df['aa']))

def make_sspair_stats(df:pd.DataFrame) -> pd.DataFrame:
    """
    Funkcija prima kao argument DataFrame za odredjeni tip rastojanja i vraca po paru ss-ss statistike.
    """
    ss_pair_df = df[["ss1", "ss2", "rastojanje"]].copy()
    ss_pair_df["rastojanje"] = pd.to_numeric(ss_pair_df["rastojanje"], errors="coerce")
    ss_pair_df = ss_pair_df.dropna(subset=["rastojanje", "ss1", "ss2"])

    ss_pair_df["ss_pair"] = list(map(tuple, np.sort(ss_pair_df[["ss1", "ss2"]].values, axis=1)))

    stats_df = (
        ss_pair_df
        .groupby("ss_pair")["rastojanje"]
        .agg(count="count", mean="mean", min="min", max="max", std="std", median="median")
        .reset_index()
    )

    stats_df["ss_pair"] = stats_df["ss_pair"].apply(lambda t: f"{t[0]} - {t[1]}")
    stats_df = stats_df.sort_values("count", ascending=False)

    return stats_df


def make_segment_hist_png(segment_lengths: list, protein: str, bins: int = 30) -> bytes:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(segment_lengths, bins=bins, color="#1f77b4", edgecolor="white")
    ax.set_title(f"Raspodela dužina segmenata — {protein}")
    ax.set_xlabel("Dužina segmenta (br. rezidua)")
    ax.set_ylabel("Broj pojavljivanja")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()
