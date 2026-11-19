import os, re, json
import pandas as pd
import build_cdh as B
import content_brsh as C
N = C.N; OUT = B.OUT; AN = B.AN; f3 = B.f3; H = C.HYP; L = C.LABEL
def lab(v): return L.get(v, v)
def esm():
    J = B.Doc(cite="ay"); J.para(f"<b>Supplementary material for: {C.TITLE}</b>")
    META = pd.read_csv(os.path.join(AN, "panel_META_BRSH_std3.csv"), encoding="utf-8-sig").set_index("variable")
    SRC = {C.W: "KCHS walking ≥30 min on ≥5 days (standardised %)", C.MV: "KCHS moderate or vigorous physical activity (standardised %)", C.FAC: "Ministry of Culture, Sports and Tourism public sports facilities",
           C.SM: "KCHS current smoking (standardised %)", C.DR: "KCHS high-risk drinking among drinkers (standardised %)", C.ED: "Census share with tertiary education (interpolated)",
           C.LIV: "Basic livelihood security benefit recipients (%)", C.DEN: "Dentists per 100,000 residents", C.MF: "Census on establishments: manufacturing workers / all workers (%)",
           "urban_population_share": "Urban-zone population share (%)", "clinics_per_100k": "Clinics per 100,000", "sex_ratio": "Males per 100 females", "fiscal_independence": "Local fiscal independence (%)",
           "pop80_share": "Population aged ≥80 years (%)", "fires_per_10000": "Fires per 10,000 residents", "multicultural_marriage_share": "Multicultural marriages (%)", "rear_seatbelt_wearing_rate": "KCHS rear seat-belt use (standardised %)"}
    rows = []
    for v in H + B.FORCED + ["fires_per_10000", "multicultural_marriage_share", "rear_seatbelt_wearing_rate"]:
        m = META.loc[v]; imp = int(m.n_imputed_interp + m.n_imputed_sgg_mean + m.n_imputed_year_mean)
        rows.append([lab(v), SRC[v], {"none": "None", "log": "Log", "per100k": "Per 100,000", "pct_of_pop": "% of population", "ratio_pop": "% of population"}.get(str(m["transform"]), str(m["transform"])), f"{m.obs_year_min}–{m.obs_year_max}", str(imp), "" if pd.isna(m.vif_final) else f"{m.vif_final:.2f}"])
    J.new_page("landscape")
    J.table_fit("Table S1. Covariates: definitions, transformation, imputed district-years and variance inflation factors", ["Covariate", "Definition and source", "Transformation", "Years", "Imputed cells", "VIF"], rows,
                "Manufacturing employment (2017–2024) was carried back to 2014–2016 (687 cells) and Sejong was assigned the national median; these cells are counted in the source data file, not as imputed here. KCHS, Korea Community Health Survey; VIF, variance inflation factor.", fs_max=8)
    # S2 민감도
    cols = [("Type III", lambda v: N["kh"].get(f"KH_III|{v}")), ("Lag 1 year", lambda v: N["lag"].get(f"1|{v}")), ("Lag 2 years", lambda v: N["lag"].get(f"2|{v}")),
            ("Binomial", lambda v: N["lik"].get(f"binom|{v}")), ("Beta", lambda v: N["lik"].get(f"beta|{v}")), ("Unweighted Gaussian", lambda v: N["lik"].get(f"gfree|{v}")),
            ("Crude prevalence", lambda v: N["lik"].get(f"crude_rate_FH_plus_old65|{v}")), ("Queen graph", lambda v: N["graph"].get(f"queen_graph|{v}")), ("Alternative prior", lambda v: N["graph"].get(f"pc_prior_alt|{v}"))]
    rows = [[lab(v)] + [f3(f(v)) if f(v) else "–" for _, f in cols] for v in H]
    J.new_page("landscape")
    J.table_fit("Table S2. Sensitivity analyses I: space-time interaction, lags, likelihoods, graph and prior", ["Covariate"] + [c for c, _ in cols], rows, "Odds ratios (95% credible intervals) per standard deviation. Crude prevalence model additionally adjusted for the share aged 65 years or older.", fs_max=8)
    cols = [("2014–2019", lambda v: N["win"].get(f"W1_2014-19|{v}")), ("2020–2024", lambda v: N["win"].get(f"W2_2020-24|{v}")), ("Excluding 2020–2021", lambda v: N.get("s19", {}).get(f"drop_2020_2021|{v}")),
            ("Restricted spatial regression", lambda v: N.get("rsr", {}).get(v)), ("County and density adjusted", lambda v: N.get("rural_OR", {}).get(v)), ("Influential cells excluded", lambda v: N.get("infl", {}).get(v)),
            ("Leave-one-province-out credible, n/17", None)]
    rows = [[lab(v)] + [(f3(f(v)) if f(v) else "–") if f else str(N["loo"][v]) for _, f in cols] for v in H]
    J.new_page("landscape")
    J.table_fit("Table S3. Sensitivity analyses II: periods, spatial confounding, urbanisation, influential observations and provinces", ["Covariate"] + [c for c, _ in cols], rows,
                "Odds ratios (95% credible intervals) per standard deviation. Influential cells: district-years with a standardised residual above 3 in absolute value were to be excluded; none met this criterion, so the estimates equal those of the main model.", fs_max=8)
    mv = N["mv"]; rows = [[lab(v), str(mv[v]["n"]), f"{100*mv[v]['cred']:.0f}", f"{100*mv[v]['sign']:.0f}", f"{mv[v]['med']:.3f}", f"{mv[v]['mn']:.3f}–{mv[v]['mx']:.3f}", "Yes" if v in C.ROB else "No"] for v in H]
    rows += [["__H__", "Negative controls (added to the main model)"]] + [[lab(k), "", "", "", f3(d), "", ""] for k, d in N["neg"].items()]
    J.new_page("portrait")
    J.table_fit("Table S4. Random covariate subsets (200 specifications) and negative controls", ["Covariate", "Specifications, n", "Credible, %", "Sign agreement, %", "Median OR", "OR range", "Robust"], rows,
                "Robust: credible in ≥75% of specifications with ≥95% sign agreement and credible in Knorr-Held type II and IV models. For negative controls the odds ratio (95% credible interval) is shown in the median column.", fs_max=9)
    my = pd.read_csv(os.path.join(AN, "results_std3", "S6_moran_year.csv"))
    rows = [[y, f"{g['mean']:.1f}", f"{g['P10']:.1f}", f"{g['P90']:.1f}", f"{g['diff']:.1f}", f"{g['ratio']:.2f}", f"{my.set_index('year').loc[int(y), 'moran_I']:.2f}"] for y, g in N["gap"].items()]
    J.new_page("portrait")
    J.table_fit("Table S5. District prevalence of toothbrushing after lunch by year", ["Year", "Mean, %", "10th percentile", "90th percentile", "Difference", "Ratio", "Moran's I"], rows,
                "Age-sex-standardised prevalence across 229 districts. Moran's I of logit prevalence with first-order contiguity weights.", fs_max=10)
    sc = pd.read_csv(os.path.join(AN, "results_std3", "S2_setcompare.csv"))
    NM = {"REF_forced_only": "Controls only", "S1_narrative_full": "All nine covariates (main model)", "S2_activity": "Physical activity (3)", "S3_risk_behaviour": "Smoking and drinking (2)", "S4_social_work": "Social and work conditions with dentists (4)", "S5_core": "Core set (5)", "S1b_obesity_sensitivity": "Main model plus obesity"}
    rows = [[NM[r.set], str(r.n_hyp), f"{r.WAIC:.1f}", f"{r.phi:.2f}", f"{r.sd_area:.3f}"] for r in sc.itertuples()]
    J.table_fit("Table S6. Covariate set comparison", ["Covariate set", "Covariates, n", "WAIC", "Spatial mixing (φ)", "Spatial SD"], rows, "WAIC, widely applicable information criterion; lower is better.", fs_max=10)
    f = os.path.join(OUT, "Supplementary_Material.docx"); J.save(f); B.pdf(f); return f
STROBE = [("Title and abstract", "1", "Design stated in title (district panel) and abstract (ecological panel study)"), ("Background/rationale", "2", "Introduction, paragraphs 1–3"), ("Objectives", "3", "Introduction, final paragraph"),
          ("Study design", "4", "Method, paragraph 1"), ("Setting", "5", "Method, paragraph 1 (229 districts, 2014–2024)"), ("Participants", "6", "Method, paragraphs 1–2 (survey respondents aggregated by district)"),
          ("Variables", "7", "Method, paragraphs 2–3; Supplementary Table S1"), ("Data sources/measurement", "8", "Method, paragraphs 1–3"), ("Bias", "9", "Method, paragraph 5 (negative controls, sensitivity analyses)"),
          ("Study size", "10", "Method, paragraph 2 (2,519 district-years)"), ("Quantitative variables", "11", "Method, paragraph 3 (standardisation, transformations)"), ("Statistical methods", "12", "Method, paragraphs 4–5"),
          ("Participants", "13", "Results, paragraph 1; Table 1"), ("Descriptive data", "14", "Table 1; Supplementary Table S5"), ("Outcome data", "15", "Results, paragraph 1; Figure 3"),
          ("Main results", "16", "Results, paragraph 2; Table 2; Figure 2"), ("Other analyses", "17", "Results, paragraph 3; Supplementary Tables S2–S4, S6"), ("Key results", "18", "Discussion, paragraph 1"),
          ("Limitations", "19", "Discussion, paragraph 5"), ("Interpretation", "20", "Discussion, paragraphs 2–4"), ("Generalisability", "21", "Discussion, paragraph 5"), ("Funding", "22", "Title page, Funding statement")]
def strobe():
    J = B.Doc(cite="ay"); J.para("<b>STROBE checklist (cross-sectional and ecological panel design)</b>")
    J.table_fit("STROBE items and location in the manuscript", ["Section/item", "No.", "Location"], [list(r) for r in STROBE], "", fs_max=10)
    f = os.path.join(OUT, "STROBE_checklist.docx"); J.save(f); B.pdf(f); return f
COVER = [
 "Dear Dr Sabbah,",
 f"Please consider the manuscript entitled “{C.TITLE}” for publication as an Original Research Article in Community Dental Health.",
 "Toothbrushing after lunch is a national oral health indicator in Korea because it depends on schools and workplaces rather than on home routines. This study follows the indicator in all 229 Korean districts from 2014 to 2024 with a small-area model that accounts for survey sampling error, and asks whether districts where adults brush after lunch are also those with other favourable health behaviours.",
 "Districts with more walking, less smoking and high-risk drinking and more tertiary education had more toothbrushing after lunch, and these associations withstood space-time interaction models, random covariate subsets and leaving out each province. Dentist supply, sports facilities, manufacturing employment and deprivation were not associated. A negative-control behaviour, rear seat-belt use, was also associated, which the manuscript reports as evidence of a general district pattern of protective behaviour. The findings speak to the journal's dental public health readership by suggesting that oral hygiene promotion could be combined with general health promotion in districts where several behaviours lag.",
 "The manuscript is original, has not been published and is not under consideration elsewhere, and will not be published elsewhere in the same form if accepted. The study used only published aggregate statistics, and the Institutional Review Board waived review. The author declares no conflicts of interest and received no funding. The analysis code is publicly archived. Generative artificial intelligence tools were used only for translation and English-language editing.",
 "Sincerely,", "Seongdae Kim, PhD", B.AFFIL, f"Email: {B.EMAIL}",
]
def cover():
    J = B.Doc(cite="ay")
    for i, p in enumerate(COVER): J.para(p, before=0 if i else 0)
    f = os.path.join(OUT, "Cover_letter.docx"); J.save(f); B.pdf(f)
    open(os.path.join(OUT, "Cover_letter_paste.txt"), "w").write("\n\n".join(COVER)); return f
if __name__ == "__main__":
    import sanitize_docx
    for fn, an in ((esm, True), (strobe, True), (cover, False)):
        f = fn(); sanitize_docx.sanitize(f, anon=an); print("[esm]", f)
