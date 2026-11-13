#!/usr/bin/env python3
"""Serum targeted metabolomics in a chronic kidney disease (CKD) exercise study: Normal vs CKD vs CKD + activity wheel (n = 10 each).
Welch t-tests on log2 concentrations, Benjamini-Hochberg FDR across metabolites per contrast,
Mann-Whitney as a robustness check, metabolite-class summaries, and a 'rescue' analysis asking
whether wheel running moves disease-altered metabolites back toward normal.

The data are unpublished and are read from data/private/ (not distributed). Public outputs are the
summary figure and aggregate statistics; per-metabolite tables are written to results/private/."""
import json, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
PRIV, RPRIV, FIG = ROOT / "data" / "private", ROOT / "results" / "private", ROOT / "results" / "figures"
RPRIV.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
if not (PRIV / "Kstudy_serum.csv").exists():
    sys.exit("Missing data/private/Kstudy_serum.csv. The data are unpublished and not distributed; see the README.")
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight"})
k = pd.read_csv(PRIV / "Kstudy_serum.csv")
grp = k.Group.map({"Con_no": "Normal", "CKD_no": "CKD", "CDK_AW": "CKD+wheel"})
X = k.iloc[:, 7:].apply(pd.to_numeric); X.index = k["Sample ID"]
const = X.columns[X.nunique() <= 1]; X = X.drop(columns=const)
L = np.log2(X.clip(lower=X[X > 0].min().min() / 2))
def cls(m):
    if re.match(r"^C\d", m): return "acylcarnitines"
    if m.startswith("lysoPC"): return "lysophosphatidylcholines"
    if m.startswith("PC"): return "phosphatidylcholines"
    if m.startswith("SM"): return "sphingomyelins"
    if m == "H1": return "hexose"
    if m in "Ala Arg Asn Asp Cit Gln Glu Gly His Ile Leu Lys Met Orn Phe Pro Ser Thr Trp Tyr Val".split(): return "amino acids"
    return "biogenic amines"
C = pd.Series({m: cls(m) for m in L.columns})
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out
def contrast(a, b):
    A, B = L[(grp == a).values], L[(grp == b).values]
    t, p = stats.ttest_ind(A, B, equal_var=False)
    mw = [stats.mannwhitneyu(A[c], B[c]).pvalue for c in L.columns]
    r = pd.DataFrame({"class": C, "log2fc": A.mean() - B.mean(), "p_welch": p, "fdr": bh(p), "p_mannwhitney": mw})
    return r
R = {"CKD vs Normal": contrast("CKD", "Normal"), "CKD+wheel vs CKD": contrast("CKD+wheel", "CKD"), "CKD+wheel vs Normal": contrast("CKD+wheel", "Normal")}
S = {"samples": int(len(X)), "per_group": grp.value_counts().to_dict(), "metabolites_analysed": int(L.shape[1]), "constant_metabolites_dropped": list(const),
     "class_sizes": C.value_counts().to_dict(), "contrasts": {}}
for n, r in R.items():
    sig = r.fdr < 0.05
    S["contrasts"][n] = {"fdr_lt_0.05": int(sig.sum()), "up": int((sig & (r.log2fc > 0)).sum()), "down": int((sig & (r.log2fc < 0)).sum()),
                         "also_mannwhitney_p_lt_0.05": int((sig & (r.p_mannwhitney < 0.05)).sum()), "nominal_p_lt_0.05": int((r.p_welch < 0.05).sum()),
                         "by_class": {c: {"n": int((C == c).sum()), "significant": int((sig & (C == c)).sum()), "mean_log2fc": round(float(r.log2fc[C == c].mean()), 2)} for c in C.unique()},
                         "top10": [{"metabolite": m, "log2fc": round(float(x.log2fc), 2), "fdr": float(f"{x.fdr:.2g}")} for m, x in r.sort_values("p_welch").head(10).iterrows()]}
d, e = R["CKD vs Normal"], R["CKD+wheel vs CKD"]
dis = d.fdr < 0.05
S["rescue"] = {"disease_altered_metabolites": int(dis.sum()),
               "correlation_disease_effect_vs_exercise_effect_all": round(float(d.log2fc.corr(e.log2fc)), 2),
               "correlation_among_disease_altered": round(float(d.log2fc[dis].corr(e.log2fc[dis])), 2) if dis.sum() > 2 else None,
               "disease_altered_moved_toward_normal": int((np.sign(e.log2fc[dis]) == -np.sign(d.log2fc[dis])).sum()),
               "median_fraction_of_disease_effect_reversed": round(float((-e.log2fc[dis] / d.log2fc[dis]).median()), 2) if dis.sum() else None,
               "still_different_from_normal_after_wheel_fdr": int(((R["CKD+wheel vs Normal"].fdr < 0.05) & dis).sum()),
               "disease_altered_and_significantly_changed_by_wheel": int((dis & (e.fdr < 0.05)).sum())}
Z = (L - L.mean()) / L.std(); u, s, _ = np.linalg.svd(Z.values, full_matrices=False); pc = u[:, :2] * s[:2]; ve = (s ** 2 / (s ** 2).sum())[:2]
S["pca_variance_explained"] = [round(float(v), 3) for v in ve]
col = {"Normal": "#8c8c8c", "CKD": "#D55E00", "CKD+wheel": "#0072B2"}
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
for g in col: ax[0].scatter(pc[(grp == g).values, 0], pc[(grp == g).values, 1], s=28, color=col[g], label=g, lw=0)
ax[0].set_xlabel(f"PC1 ({ve[0]:.0%})"); ax[0].set_ylabel(f"PC2 ({ve[1]:.0%})"); ax[0].legend(frameon=False, fontsize=8); ax[0].set_title("Serum metabolome by group", loc="left")
ax[1].scatter(d.log2fc[~dis], -np.log10(d.p_welch[~dis]), s=8, color="#bdbdbd", lw=0); ax[1].scatter(d.log2fc[dis], -np.log10(d.p_welch[dis]), s=12, color="#D55E00", lw=0)
ax[1].set_xlabel("log2 fold change, CKD vs Normal"); ax[1].set_ylabel("-log10 p"); ax[1].set_title(f"Disease effect: {int(dis.sum())} of {len(d)} metabolites at FDR < 0.05", loc="left")
ax[2].scatter(d.log2fc[~dis], e.log2fc[~dis], s=8, color="#bdbdbd", lw=0); ax[2].scatter(d.log2fc[dis], e.log2fc[dis], s=12, color="#0072B2", lw=0, label="altered by CKD")
lim = max(d.log2fc.abs().max(), e.log2fc.abs().max()) * 1.05; ax[2].plot([-lim, lim], [lim, -lim], color="#8c8c8c", ls="--", lw=.8, label="complete reversal")
ax[2].axhline(0, color="#8c8c8c", lw=.6); ax[2].axvline(0, color="#8c8c8c", lw=.6); ax[2].legend(frameon=False, fontsize=7)
ax[2].set_xlabel("disease effect (CKD vs Normal, log2)"); ax[2].set_ylabel("exercise effect (CKD+wheel vs CKD, log2)"); ax[2].set_title("Does wheel running reverse the disease signature?", loc="left")
fig.tight_layout(); fig.savefig(FIG / "fig1_serum_metabolome.png"); plt.close(fig)
with pd.ExcelWriter(RPRIV / "per_metabolite_results.xlsx") as xw:
    for n, r in R.items(): r.sort_values("p_welch").round(5).to_excel(xw, sheet_name=n.replace("+", " plus ")[:31])
json.dump(S, open(RPRIV / "summary_with_top_metabolites.json", "w"), indent=2)
public = json.loads(json.dumps(S)); public.pop("constant_metabolites_dropped")
public["constant_metabolites_dropped"] = len(S["constant_metabolites_dropped"])
for c in public["contrasts"].values():
    c.pop("top10")                     # metabolite-level results stay private
json.dump(public, open(ROOT / "results" / "summary.json", "w"), indent=2)
out = {k_: v for k_, v in S.items() if k_ != "contrasts"}; print(json.dumps(out))
for n, v in S["contrasts"].items():
    print(n, {a: b for a, b in v.items() if a not in ("by_class", "top10")}); print("  ", {c: (x["significant"], x["n"], x["mean_log2fc"]) for c, x in v["by_class"].items()}); print("  ", [(t["metabolite"], t["log2fc"], t["fdr"]) for t in v["top10"]])
