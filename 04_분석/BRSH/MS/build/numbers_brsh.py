#!/usr/bin/env python3
import os, json
import pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); AN = os.path.normpath(os.path.join(HERE, "..", "..")); R = os.path.join(AN, "results_std3")
ROOT = os.path.normpath(os.path.join(AN, "..", ".."))
rd = lambda f: pd.read_csv(os.path.join(R, f))
H = ["walking_practice_rate", "moderate_pa_rate", "public_sports_facilities_per_100k", "current_smoking_rate", "high_risk_drinking_rate_drinkers",
     "tertiary_educated_share", "basic_livelihood_rate", "dentists_per_100k", "ext__manufacturing_share"]
F = ["urban_population_share", "clinics_per_100k", "sex_ratio", "fiscal_independence", "pop80_share"]
d = lambda r: dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), cred=bool(r.cred))
N = {}
m = rd("S3_main_OR.csv"); N["main"] = {r.cov: d(r) for r in m.itertuples()}
fit = rd("S3_main_fit.csv").iloc[0]; N["fit"] = {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v)) for k, v in fit.items()}
s0 = rd("S0_screen.csv").iloc[0]; N["moran"] = float(s0.moran_I_pooled_logitSR); N["persist"] = float(s0.persistence_rank_rho); N["p9010"] = float(s0.P90_P10_SR)
N["phi_forced"] = float(s0.phi_forced_only)
my = rd("S6_moran_year.csv"); N["moran_ymin"] = float(my.moran_I.min()); N["moran_ymax"] = float(my.moran_I.max()); N["kendallW"] = float(rd("S6_kendallW.csv").kendall_W[0])
hp = rd("S3_main_hyper.csv").set_index("hyperpar"); N["phi"] = float(hp.loc["Phi for area", "mean"]); N["phi_lo"] = float(hp.loc["Phi for area", "X0.025quant"]); N["phi_hi"] = float(hp.loc["Phi for area", "X0.975quant"])
vp = rd("S3_varpart.csv").set_index("component"); N["vp"] = {k: float(v) for k, v in vp.share.items()}
nt = vp.drop("temporal").variance; N["vp_notime"] = {k: float(100 * v / nt.sum()) for k, v in nt.items()}
sp = rd("S3_spatial_OR.csv"); N["n_hot"] = int((sp["class"] == "hotspot").sum()); N["n_cold"] = int((sp["class"] == "coldspot").sum())
N["hot_sido"] = sp[sp["class"] == "hotspot"].sido.value_counts().to_dict(); N["cold_sido"] = sp[sp["class"] == "coldspot"].sido.value_counts().to_dict()
sp0 = rd("S0_forced_only_spatialOR.csv"); N["hot0_sido"] = sp0[sp0["class"] == "hotspot"].sido.value_counts().to_dict(); N["cold0_sido"] = sp0[sp0["class"] == "coldspot"].sido.value_counts().to_dict()
N["n_hot0"] = int((sp0["class"] == "hotspot").sum()); N["n_cold0"] = int((sp0["class"] == "coldspot").sum())
N["vs_forced"] = rd("S3_vs_forced_only.csv").iloc[0].to_dict()
N["neg"] = {r.cov: d(r) for r in rd("S3_negcontrol.csv").itertuples()}
k = rd("S8_KH.csv"); N["kh"] = {f"{r.model}|{r.cov}": d(r) for r in k.itertuples()}; N["kh_waic"] = k.groupby("model").WAIC.first().to_dict(); N["kh_waic"]["main"] = float(fit.WAIC)
b = rd("S9_between_within.csv"); N["bw"] = {f"{r.cov}|{r.component}": d(r) for r in b.itertuples()}
w = rd("S4_window_OR.csv"); N["win"] = {f"{r.window}|{r.cov}": d(r) for r in w.itertuples()}
ia = rd("S4_interaction.csv"); N["win_ix"] = {f"{r.window}|{r.cov}": d(r) for r in ia.itertuples()}
cv = rd("S5_covid.csv"); N["covid"] = {r.cov: d(r) for r in cv.itertuples()}
lg = rd("S10_lag.csv"); N["lag"] = {f"{r.lag}|{r.cov}": d(r) for r in lg.itertuples()}
def opt(f, key, fn):
    p = os.path.join(R, f)
    if os.path.exists(p): N[key] = fn(pd.read_csv(p))
opt("S7_multiverse_summary.csv", "mv", lambda x: {r.cov: dict(n=int(r.n_specs), cred=float(r.cred_frac), sign=float(r.sign_agree), med=float(r.OR_median), mn=float(r.OR_min), mx=float(r.OR_max), robust=bool(r.robust)) for r in x.itertuples()})
if "mv" not in N and os.environ.get("MV_PARTIAL"):   # 시험 빌드 전용(최종 아님)
    import glob as _g
    mvp = pd.concat([pd.read_csv(f) for f in _g.glob(os.path.join(R, "S7_multiverse_results_part_*.csv"))]); sg = {r.cov: np.sign(np.log(r.OR)) for r in m.itertuples()}
    N["mv"] = {v: dict(n=int((mvp["cov"] == v).sum()), cred=float(mvp[mvp["cov"] == v].cred.mean()), sign=float((np.sign(np.log(mvp[mvp["cov"] == v].OR)) == sg[v]).mean()), med=float(mvp[mvp["cov"] == v].OR.median()), mn=float(mvp[mvp["cov"] == v].OR.min()), mx=float(mvp[mvp["cov"] == v].OR.max()), robust=bool(mvp[mvp["cov"] == v].cred.mean() >= .75 and (np.sign(np.log(mvp[mvp["cov"] == v].OR)) == sg[v]).mean() >= .95)) for v in H}
    N["mv_partial"] = True
opt("S11_leave_one_sido.csv", "loo", lambda x: {c: int(x[x["cov"] == c].cred.sum()) for c in H}); opt("S11_leave_one_sido.csv", "loo_n", lambda x: int(x.left_out.nunique()))
opt("S12_likelihood.csv", "lik", lambda x: {f"{r.likelihood}|{r.cov}": d(r) for r in x.itertuples()})
opt("S13_graph_prior.csv", "graph", lambda x: {f"{r.variant}|{r.cov}": d(r) for r in x.itertuples()})
opt("S14_rural_OR.csv", "rural_OR", lambda x: {r.cov: d(r) for r in x.itertuples()}); opt("S14_rural_summary.csv", "rural", lambda x: x.iloc[0].to_dict())
opt("S16_RSR.csv", "rsr", lambda x: {r.cov: d(r) for r in x.itertuples()} if "OR" in x else {})
opt("S17_region_slopes.csv", "region", lambda x: {f"{r.covariate}|{r.region}": d(r) for r in x.itertuples()} if "OR" in x else {})
opt("S18_influence.csv", "infl", lambda x: {r.cov: d(r) for r in x.itertuples()}); opt("S18_influence.csv", "infl_variant", lambda x: str(x.variant[0]))
opt("S19_breaks.csv", "s19", lambda x: {f"{r.variant}|{r.cov}": d(r) for r in x.itertuples()})
opt("S20_gap_trend.csv", "gap", lambda x: {str(int(r.year)): dict(mean=float(r.mean_SR), P90=float(r.P90), P10=float(r.P10), ratio=float(r.P90_P10), diff=float(r.P90_minus_P10), SD=float(r.SD_logit)) for r in x.itertuples()})
opt("S20_slope_summary.csv", "slope", lambda x: x.iloc[0].to_dict())
FO = pd.read_csv(os.path.join(ROOT, "02_패널", "BRSH_fh_outcomes.csv"), encoding="utf-8-sig")
N["sr"] = {str(y): round(float(v), 1) for y, v in FO.groupby("year").SR.mean().items()}; N["sr_mean_all"] = round(float(FO.SR.mean()), 1)
N["n_cells"] = int(len(FO)); N["n_resp_total"] = int(FO.N.sum()); N["n_resp_median"] = int(FO.N.median())
pv = FO.pivot(index="fbd_region", columns="year", values="SR")
N["rank_rho"] = {f"{a}-{a+1}": round(float(pv[a].corr(pv[a + 1], method="spearman")), 2) for a in range(2014, 2024)}
N["sido_raw"] = FO.assign(s=FO.fbd_region.str[:3]).groupby("s").SR.mean().round(1).to_dict()
if "gap" not in N:   # S20 이전 시험 빌드: 관측값에서 같은 정의로 계산
    lg = np.log(FO.SR / (100 - FO.SR))
    N["gap"] = {str(y): dict(mean=float(g.SR.mean()), P90=float(g.SR.quantile(.9)), P10=float(g.SR.quantile(.1)), ratio=float(g.SR.quantile(.9) / g.SR.quantile(.1)), diff=float(g.SR.quantile(.9) - g.SR.quantile(.1)), SD=float(lg[g.index].std())) for y, g in FO.groupby("year")}
mv = N.get("mv", {}); kh = N["kh"]
N["robust"] = [v for v in H if mv.get(v, {}).get("robust") and all(kh.get(f"{t}|{v}", {}).get("cred") for t in ("KH_II", "KH_IV"))]
json.dump(N, open(os.path.join(HERE, "numbers_brsh.json"), "w"), ensure_ascii=False, indent=1, default=float)
print("robust:", N["robust"]); print({k: (len(v) if isinstance(v, dict) else v) for k, v in N.items() if k in ("mv", "loo", "lik", "rsr", "s19", "gap", "rural")})
