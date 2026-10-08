"""Catalogue sources of a run (R16e, C1): the figures of the primary family and the QC figures, built from the run's own outputs.

The sources are written as TSV tables with the columns that the catalogue builders read. Their values are copied from the run:
the analysis matrix and its observed mask (preprocessing/primary), the observations (groups), the PCA scores and variance
(preprocessing/qc), and the DE table of each contrast of the primary family (dea/zero_null.tsv). The only derived values are the
group means of the MA plot (means of the analysis matrix over the observed cells), the list of DEPs (q at or below the family
cutoff), and the top DEPs by q for the dot plot. No statistic is recomputed.

``build_jobs`` returns the jobs (figure kind, output ID, source path, options) and the figures that were skipped with a reason.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

PRIMARY = "primary"
TOP_DEPS_DOTPLOT = 6
MAX_VOLCANO_LABELS = 10


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, header: list[str], rows: list[list[Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    for row in rows:
        cells = []
        for value in row:
            if value is None:
                cells.append("NA")
            elif isinstance(value, bool):
                cells.append("TRUE" if value else "FALSE")
            elif isinstance(value, float):
                cells.append(repr(value))
            else:
                cells.append(str(value))
        lines.append("\t".join(cells))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _num(text: str | None) -> float | None:
    if text is None or text in ("", "NA", "NaN"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _matrix(run: Path) -> tuple[list[str], dict[str, dict[str, float | None]], dict[str, dict[str, bool]], dict[str, str]]:
    """Samples, values and observed mask of the analysis matrix, and the group of each sample (preprocessing/primary)."""
    primary = run / "preprocessing" / "primary"
    matrix = _read(primary / "matrix.tsv")
    mask = {r["feature_id"]: r for r in _read(primary / "observed_mask.tsv")}
    samples = [c for c in matrix[0] if c != "feature_id"]
    values = {r["feature_id"]: {s: _num(r[s]) for s in samples} for r in matrix}
    observed = {f: {s: mask[f][s] == "true" for s in samples} for f in values}
    groups = {r["observation_id"]: r["group"] for r in _read(primary / "observations.tsv")}
    return samples, values, observed, groups


def _contrast_plan(run: Path) -> tuple[list[dict[str, Any]], dict[str, dict[str, str]], dict[str, float]]:
    """The primary-family contrasts (in the order of the configuration), their numerator and denominator groups and labels,
    the DE rows of each contrast, and the q cutoff of each family."""
    config = json.loads((run / "config.resolved.json").read_text(encoding="utf-8"))
    families = {r["family_id"]: r for r in _read(run / "dea" / "families.tsv")}
    de_rows = [r for r in _read(run / "dea" / "zero_null.tsv") if r["eligibility"] == "tested"]
    cutoffs = {fid: float(r["q_cutoff"]) for fid, r in families.items()}
    plan = []
    for contrast in config.get("contrasts", []):
        cid = contrast["id"]
        rows = [r for r in de_rows if r["contrast_id"] == cid]
        if not rows:
            continue
        family = rows[0]["family_id"]
        if families.get(family, {}).get("role") != PRIMARY:
            continue
        weights = contrast.get("weights", {})
        numerator = [k.split(".", 1)[1] for k, w in weights.items() if w > 0]
        denominator = [k.split(".", 1)[1] for k, w in weights.items() if w < 0]
        if len(numerator) != 1 or len(denominator) != 1:
            continue
        plan.append({"id": cid, "label": contrast.get("label", cid), "family": family, "numerator": numerator[0],
                     "denominator": denominator[0], "rows": rows, "cutoff": cutoffs[family]})
    return plan, {}, cutoffs


def _qc_jobs(run: Path, work: Path) -> list[dict[str, Any]]:
    samples, values, observed, groups = _matrix(run)
    features = list(values)
    long_rows = [[s, f, values[f][s], observed[f][s]] for s in samples for f in features]
    jobs = []
    for figure in ("qc_intensity_distributions", "qc_missingness_heatmap", "qc_sample_correlation_heatmap"):
        if figure == "qc_intensity_distributions":
            header, rows = ["sample", "feature", "value", "observed"], long_rows
        elif figure == "qc_missingness_heatmap":
            header, rows = ["feature", "sample", "observed", "value"], [[f, s, observed[f][s], values[f][s]] for f in features for s in samples]
        else:
            header, rows = ["sample", "feature", "observed", "value"], [[s, f, observed[f][s], values[f][s]] for s in samples for f in features]
        jobs.append({"figure_id": figure, "output_id": figure, "source": _write(work / f"{figure}.tsv", header, rows), "options": {
            "title": {"qc_intensity_distributions": "Intensity distributions (observed cells)", "qc_missingness_heatmap": "Missingness (observed or missing)",
                      "qc_sample_correlation_heatmap": "Sample correlations (observed cells)"}[figure],
            "caption": "analysis matrix of the run; missing cells are not plotted as values"}})
    scores = {r["observation_id"]: r for r in _read(run / "preprocessing" / "qc" / "pca_scores.tsv")}
    variance = {r["component"]: float(r["variance_explained"]) * 100.0 for r in _read(run / "preprocessing" / "qc" / "pca_variance.tsv")}
    rows = [[s, groups[s], float(scores[s]["PC1"]), float(scores[s]["PC2"]), variance["PC1"], variance["PC2"]] for s in samples]
    jobs.append({"figure_id": "qc_pca_ellipses", "output_id": "qc_pca_ellipses", "source": _write(work / "qc_pca_ellipses.tsv",
                 ["sample", "group", "PC1", "PC2", "PC1_var_pct", "PC2_var_pct"], rows),
                 "options": {"title": "PCA with 95% ellipse per group", "caption": "PCA of the display matrix of the run; percent variance per axis"}})
    return jobs


def _permanova_jobs(run: Path, work: Path) -> list[dict[str, Any]]:
    """PCoA with centroids and PERMDISP, per feature set of the PERMANOVA stage (permanova/*.tsv, written by the producer)."""
    base = run / "permanova"
    if not (base / "ordination.tsv").is_file():
        return []
    ordination = _read(base / "ordination.tsv")
    axes = _read(base / "ordination_axes.tsv")
    tests = _read(base / "tests.tsv")
    jobs = []
    for set_id in dict.fromkeys(r["feature_set_id"] for r in ordination):
        members = [r for r in ordination if r["feature_set_id"] == set_id]
        pct = {r["axis"]: float(r["percent_positive_eigenvalue"]) for r in axes if r["feature_set_id"] == set_id}
        overall = next(r for r in tests if r["feature_set_id"] == set_id and r["analysis"] == "global")
        caption = f"PERMANOVA P = {overall['p_display']}; B = {overall['nperm']}"
        jobs.append({"figure_id": "perm_pcoa_centroids", "output_id": f"perm_pcoa_centroids__{set_id}",
                     "source": _write(work / f"pcoa_{set_id}.tsv", ["observation_id", "group", "axis1", "axis2", "percent_axis1", "percent_axis2"],
                                      [[r["observation_id"], r["group"], float(r["axis1"]), float(r["axis2"]), pct["axis1"], pct["axis2"]] for r in members]),
                     "options": {"title": f"Principal coordinates, feature set {set_id}", "caption": caption}})
        jobs.append({"figure_id": "perm_permdisp", "output_id": f"perm_permdisp__{set_id}",
                     "source": _write(work / f"permdisp_{set_id}.tsv", ["observation_id", "group", "distance_to_centroid", "permdisp_p", "permutations", "seed"],
                                      [[r["observation_id"], r["group"], float(r["distance_to_centroid"]), float(overall["permdisp_p"]),
                                        int(overall["nperm"]), int(overall["seed"])] for r in members]),
                     "options": {"title": f"Distances to group centroids, feature set {set_id}",
                                 "caption": f"PERMDISP over {len(members)} units; B = {overall['nperm']}, seed {overall['seed']}"}})
    return jobs


def _core_jobs(run: Path, work: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """All run-level catalogue jobs and the figures that were skipped, with their reasons."""
    jobs: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    if (run / "preprocessing" / "primary" / "matrix.tsv").is_file():
        jobs += _qc_jobs(run, work)
    else:
        skipped.append({"figure": "qc and differential", "reason": "no analysis matrix in this run"})
        return _permanova_jobs(run, work), skipped
    jobs += _permanova_jobs(run, work)
    if not (run / "dea" / "zero_null.tsv").is_file():
        skipped.append({"figure": "differential", "reason": "no DE table of a primary family in this run"})
        return jobs, skipped
    samples, values, observed, groups = _matrix(run)
    plan, _, _ = _contrast_plan(run)
    venn_sets: list[dict[str, Any]] = []
    for contrast in plan:
        cid = contrast["id"]
        num, den = contrast["numerator"], contrast["denominator"]
        cutoff, label = contrast["cutoff"], contrast["label"]
        rows = contrast["rows"]
        in_groups = [s for s in samples if groups[s] in (num, den)]
        def counts(r: dict[str, str]) -> dict[str, int]:
            return {k: int(v) for k, v in json.loads(r["n_obs_by_required_group"]).items()}
        tested = []
        for r in rows:
            q, p = _num(r["q_value"]), _num(r["p_value"])
            if q is None or p is None:
                continue
            tested.append({"feature": r["feature_id"], "effect": float(r["effect"]), "P": p, "q": q, "n": counts(r)})
        deps = sorted([t for t in tested if t["q"] <= cutoff], key=lambda t: (t["q"], t["feature"]))
        p_threshold = max((t["P"] for t in deps), default=None)
        caption = f"family {contrast['family']} (BH), q cutoff {cutoff:g}"
        # volcano: all tested features of the contrast; labels for the smallest q
        labels = [t["feature"] for t in sorted(tested, key=lambda t: (t["q"], t["feature"]))[:MAX_VOLCANO_LABELS]]
        vrows = [[t["feature"], t["effect"], t["P"], t["q"], t["n"][num], t["n"][den]] for t in tested]
        jobs.append({"figure_id": "diff_volcano_labelled", "output_id": f"diff_volcano_labelled__{cid}",
                     "source": _write(work / f"volcano_{cid}.tsv", ["feature", "effect", "P", "q", f"n_{num}", f"n_{den}"], vrows),
                     "options": {"labels": labels, "q_cutoff": cutoff, "p_threshold": p_threshold, "groups": [num, den],
                                 "title": label, "caption": caption}})
        # MA plot: group means of the analysis matrix (observed cells) against the DE values
        mrows = []
        for t in tested:
            f = t["feature"]
            mean_num = _mean([values[f][s] for s in samples if groups[s] == num and observed[f][s] and values[f][s] is not None])
            mean_den = _mean([values[f][s] for s in samples if groups[s] == den and observed[f][s] and values[f][s] is not None])
            if mean_num is None or mean_den is None:
                continue
            mrows.append([cid, num, den, f, mean_num, mean_den, t["P"], t["q"], t["n"][num], t["n"][den]])
        if mrows:
            jobs.append({"figure_id": "diff_ma_plot", "output_id": f"diff_ma_plot__{cid}",
                         "source": _write(work / f"ma_{cid}.tsv", ["contrast", "numerator", "denominator", "feature", "mean_numerator",
                                                                   "mean_denominator", "P", "q", "n_numerator", "n_denominator"], mrows),
                         "options": {"title": label, "caption": caption}})
        else:
            skipped.append({"figure": f"diff_ma_plot__{cid}", "reason": "no feature with observed values in both groups"})
        # DEP heatmap and dot plot
        dep_features = [t["feature"] for t in deps if all(observed[t["feature"]][s] and values[t["feature"]][s] is not None for s in in_groups)]
        if len(dep_features) < 2:
            skipped.append({"figure": f"diff_dep_heatmap_clustered__{cid}", "reason": f"{len(dep_features)} DEP(s) at the family cutoff; a heatmap needs at least 2"})
        else:
            hrows = [[f, s, values[f][s], True, groups[s]] for f in dep_features for s in in_groups]
            jobs.append({"figure_id": "diff_dep_heatmap_clustered", "output_id": f"diff_dep_heatmap_clustered__{cid}",
                         "source": _write(work / f"dep_{cid}.tsv", ["feature", "sample", "value", "declared", "group"], hrows),
                         "options": {"title": label, "caption": caption + "; row z-scores of the DEPs"}})
        top = [t for t in deps if all(observed[t["feature"]][s] and values[t["feature"]][s] is not None for s in in_groups)][:TOP_DEPS_DOTPLOT]
        if not top:
            skipped.append({"figure": f"diff_protein_dotplot_brackets__{cid}", "reason": "no DEP at the family cutoff"})
        else:
            drows = [[t["feature"], groups[s], s, values[t["feature"]][s]] for t in top for s in in_groups]
            family = [{"protein": t["feature"], "comparison": f"{num} - {den}", "group_a": den, "group_b": num, "effect": t["effect"],
                       "P": t["P"], "q": t["q"], f"n_{num}": t["n"][num], f"n_{den}": t["n"][den]} for t in top]
            jobs.append({"figure_id": "diff_protein_dotplot_brackets", "output_id": f"diff_protein_dotplot_brackets__{cid}",
                         "source": _write(work / f"dotplot_{cid}.tsv", ["protein", "group", "replicate", "value"],
                                          [[p, g, s, v] for p, g, s, v in drows]),
                         "options": {"family": family, "title": label,
                                     "caption": f"top {len(top)} DEPs by q; brackets: family P"}})
        venn_sets.append({"id": cid, "label": label, "deps": {t["feature"] for t in deps}, "tested": {t["feature"] for t in tested}})
    if len(venn_sets) >= 2:
        universe = sorted(set().union(*(s["tested"] for s in venn_sets)))
        letters = [chr(ord("A") + i) for i in range(len(venn_sets))]
        vrows = [[f] + [1 if f in s["deps"] else 0 for s in venn_sets] for f in universe]
        jobs.append({"figure_id": "diff_upset_venn", "output_id": "diff_upset_venn__primary",
                     "source": _write(work / "venn_primary.tsv", ["feature"] + letters, vrows),
                     "options": {"sets": letters, "set_labels": [s["label"] for s in venn_sets], "venn": len(letters) <= 3,
                                 "title": "DEPs across the primary contrasts", "caption": "DEPs at the family cutoff"}})
    else:
        skipped.append({"figure": "diff_upset_venn__primary", "reason": "fewer than two primary contrasts"})
    return jobs, skipped


def _mean(values: list[float | None]) -> float | None:
    present = [v for v in values if v is not None]
    return sum(present) / len(present) if present else None


def _ppi_snapshot(run: Path) -> Path | None:
    """The interaction snapshot of the run's PPI resource (its edge file, named in the snapshot manifest)."""
    if not (run / "config.resolved.json").is_file():
        return None
    config = json.loads((run / "config.resolved.json").read_text(encoding="utf-8"))
    for resource in config.get("resources", []):
        if resource.get("kind") == "ppi":
            manifest = Path(resource["path"])
            if not manifest.is_file():
                return None
            for item in json.loads(manifest.read_text(encoding="utf-8")).get("files", []):
                if item.get("role") == "ppi":
                    return manifest.parent / item["name"]
    return None


def _postde_jobs(run: Path, work: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Post-DE figures, from the result tables of the modules that ran (no statistic is recomputed; the derived columns are the
    ones named in the receipt: the pooled AUC is the producer's, the permutation P is the contract rule on the published null,
    the module eigengenes are group means of the producer's values)."""
    jobs: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    biomarker = run / "post_de" / "biomarker"
    figure_sources = biomarker / "figure_sources"
    if (biomarker / "cv_performance.tsv").is_file() and (figure_sources / "roc_cv.tsv").is_file() and (figure_sources / "roc_cv_band.tsv").is_file():
        performance = _read(biomarker / "cv_performance.tsv")[0]
        claim, seed = performance["claim_label"], performance["seed"]
        nulls = _read(figure_sources / "permutation_null.tsv")
        band_draws = min(500, int(performance["threshold_bootstrap_draws"]))   # the producer's ROC band rule (min(500, bootstrap))
        counts = [performance["n_negative_units"], performance["n_positive_units"]]
        header = ["component", "fpr", "sensitivity", "band_lower", "band_upper", "auc", "auc_lower", "auc_upper", "n_control", "n_case",
                  "claim_label", "score_source", "seed", "B"]
        common = [performance["pooled_oof_auc"], "NA", "NA", counts[0], counts[1], claim, "out_of_fold", seed, band_draws]
        roc = [["curve", r["fpr"], r["tpr"], "NA", "NA"] for r in _read(figure_sources / "roc_cv.tsv")]
        band = [["band", r["fpr"], "NA", r["tpr_lower"], r["tpr_upper"]] for r in _read(figure_sources / "roc_cv_band.tsv")]
        jobs.append({"figure_id": "bio_roc_ci_band", "output_id": "bio_roc_ci_band", "source": _write(work / "bio_roc.tsv", header, [r + common for r in roc + band]),
                     "options": {"title": "ROC of the cross-validated classifier", "caption": f"claim label: {claim}; out-of-fold predictions"}})
        repeats = _read(figure_sources / "cv_auc_repeats.tsv")
        if performance["permutation_p"] in ("NA", "") or performance["permutation_B"] in ("NA", ""):
            skipped.append({"figure": "bio_cv_auc_vs_permutation", "reason": "the run has no permutation test (permutation disabled)"})
        else:
            # The displayed P is the producer's published value; the builder checks it against (k + 1) / (B + 1) and refuses a mismatch.
            perm_b, p_display = int(performance["permutation_B"]), float(performance["permutation_p"])
            cv_rows = [["cv_repeat", r["auc"], perm_b, seed, p_display, claim] for r in repeats] + \
                      [["null", n["pooled_oof_auc"], perm_b, seed, p_display, claim] for n in nulls]
            jobs.append({"figure_id": "bio_cv_auc_vs_permutation", "output_id": "bio_cv_auc_vs_permutation",
                         "source": _write(work / "bio_cv.tsv", ["series", "value", "B", "seed", "P_display", "claim_label"], cv_rows),
                         "options": {"title": "CV AUC against the permutation null", "caption": f"claim label: {claim}"}})
        single = sorted(_read(biomarker / "single_feature_auc.tsv"), key=lambda r: -float(r["auc"]))[:10]
        forest = [[r["feature_id"], r["auc"], r["bootstrap_lower"], r["bootstrap_upper"], r["n_negative_units"], r["n_positive_units"],
                   r["claim_label"]] for r in single]
        jobs.append({"figure_id": "bio_panel_forest", "output_id": "bio_panel_forest",
                     "source": _write(work / "bio_forest.tsv", ["feature", "estimate", "ci_lower", "ci_upper", "n_Control", "n_Case", "claim_label"], forest),
                     "options": {"xlab": "single-feature AUC (95% bootstrap interval)", "title": "Single-feature AUC, top 10 features",
                                 "caption": f"in-sample (optimistic); not a panel effect; claim label: {single[0]['claim_label']}; no P value"}})
    else:
        skipped.append({"figure": "biomarker figures", "reason": "the biomarker module did not run (no CV tables)"})
    pathways = run / "pathways"
    for direction in ("ora_up", "ora_down"):
        table = pathways / f"{direction}.tsv"
        if not table.is_file():
            skipped.append({"figure": f"path_enrichment_dotplot__{direction}", "reason": "no over-representation table (the pathway stage did not run)"})
            continue
        rows = [r for r in _read(table) if r["status"] == "tested"]
        if not any(int(r["k"]) >= 1 for r in rows):
            skipped.append({"figure": f"path_enrichment_dotplot__{direction}", "reason": "no gene set with a foreground hit"})
            continue
        header = ["set_id", "N", "K", "n", "k", "p_value", "q_value", "overlap_genes"]
        jobs.append({"figure_id": "path_enrichment_dotplot", "output_id": f"path_enrichment_dotplot__{direction}",
                     "source": _write(work / f"{direction}.tsv", header, [[r[c] for c in header] for r in rows]),
                     "options": {"title": f"Over-representation, {direction.replace('ora_', '')} sets", "caption": f"family {rows[0]['family_id']} (BH)"}})
    networks = run / "post_de" / "networks"
    observations = run / "preprocessing" / "primary" / "observations.tsv"
    if (networks / "modules.tsv").is_file() and (networks / "eigengenes.tsv").is_file() and observations.is_file():
        assigned = sorted({r["module"] for r in _read(networks / "modules.tsv") if r["module"] not in ("unassigned", "")})
        eigen = _read(networks / "eigengenes.tsv")
        modules = [m for m in assigned if m in eigen[0]]
        if not modules:
            skipped.append({"figure": "net_module_heatmap", "reason": "no co-abundance module was assigned"})
        else:
            group_of = {r["observation_id"]: r["group"] for r in _read(observations)}
            groups = list(dict.fromkeys(group_of[r["observation_id"]] for r in eigen))
            rows = []
            for m in modules:
                for g in groups:
                    values = [float(r[m]) for r in eigen if group_of[r["observation_id"]] == g]
                    rows.append([m, g, sum(values) / len(values), len(values), "TRUE"])
            jobs.append({"figure_id": "net_module_heatmap", "output_id": "net_module_heatmap",
                         "source": _write(work / "modules.tsv", ["module", "group", "eigengene", "n", "declared"], rows),
                         "options": {"modules": modules, "title": "Module eigengenes, group means",
                                     "caption": "module level, descriptive"}})
    else:
        skipped.append({"figure": "net_module_heatmap", "reason": "the network stage did not write its module tables"})
    edges_table = networks / "ppi_subnetwork.tsv"
    if edges_table.is_file():
        snapshot = _ppi_snapshot(run)
        if snapshot is None or not snapshot.is_file():
            skipped.append({"figure": "net_network_graph", "reason": "the run has no PPI resource to check the subnetwork edges against"})
        else:
            edges = _read(edges_table)
            for set_id in dict.fromkeys(e["set_id"] for e in edges):
                members = [e for e in edges if e["set_id"] == set_id]
                nodes = list(dict.fromkeys([m["gene_a"] for m in members] + [m["gene_b"] for m in members]))
                rows = [["node", n, "", "", set_id] for n in nodes] + [["edge", f"{m['gene_a']}-{m['gene_b']}", m["gene_a"], m["gene_b"], ""] for m in members]
                jobs.append({"figure_id": "net_network_graph", "output_id": f"net_network_graph__{set_id}",
                             "source": _write(work / f"network_{set_id}.tsv", ["kind", "id", "from", "to", "module"], rows),
                             "options": {"snapshot": snapshot, "title": f"Declared subnetwork ({set_id})",
                                         "caption": "edges as in the run's PPI resource; layout is presentation only"}})
    return jobs, skipped


def build_jobs(run: Path, work: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """All run-level catalogue jobs (QC, differential, PERMANOVA, post-DE) and the figures that were skipped, with reasons."""
    jobs, skipped = _core_jobs(run, work)
    post_jobs, post_skipped = _postde_jobs(run, work)
    return jobs + post_jobs, skipped + post_skipped

