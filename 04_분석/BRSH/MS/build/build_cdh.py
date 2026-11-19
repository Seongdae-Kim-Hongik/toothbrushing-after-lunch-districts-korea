#  구조 초록 ≤250(Objective/Basic research design/Participants/Main outcome measures/Results/Conclusions, 굵은 이탤릭+쌍점), 키워드 ≥5,
#  본문 Introduction·Method·Results·Discussion·Conclusion, 소제목·각주 최소, 10pt TNR 더블 A4 여백 3cm, 표는 원고 끝('Table n near here'),
#  Sage Harvard 참고문헌(원칙 ≤20), 그림 300dpi 이상 TIFF, EQUATOR 체크리스트 업로드.
import os, re, json, subprocess, sys
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from jkms_engine import JKMS
import content_brsh as C
import pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); MS = os.path.dirname(HERE); OUT = os.path.join(MS, "out"); os.makedirs(OUT, exist_ok=True)
AN = os.path.normpath(os.path.join(MS, "..")); ROOT = os.path.normpath(os.path.join(AN, "..", ".."))
N = C.N; REFS = json.load(open(os.path.join(HERE, "refs_brsh.json"), encoding="utf-8"))
rdf = lambda f, d: open(os.path.join(OUT, f)).read().strip() if os.path.exists(os.path.join(OUT, f)) else d
ZENODO = rdf("zenodo_doi.txt", "10.5281/zenodo.PENDING"); GITHUB = rdf("github_url.txt", "https://github.com/PENDING")
ORCID = "0000-0002-1525-0091"
AFFIL = "Department of Biohealth Convergence Open Sharing System, Hongik University, Seoul, Republic of Korea"
ADDR = "Department of Biohealth Convergence Open Sharing System, Hongik University, 94 Wausan-ro, Mapo-gu, Seoul 04066, Republic of Korea"
EMAIL = "kimsd@hongik.ac.kr"
JFULL = {"lee2023factors": "BMC Oral Health", "son2026health": "Gerodontology", "son2026toothbrushing": "BMC Oral Health", "chang2021better": "Journal of Periodontology",
         "sanders2005clustering": "Community Dental Health", "jordao2018clustering": "Community Dentistry and Oral Epidemiology", "raison2021toothbrushing": "Community Dentistry and Oral Epidemiology",
         "lertpimonchai2017association": "International Dental Journal", "bakri2023oral": "Health Promotion International", "peres2019oral": "The Lancet",
         "kang2015korea": "Osong Public Health and Research Perspectives", "vonelm2007strengthening": "The Lancet", "fay1979estimates": "Journal of the American Statistical Association",
         "riebler2016intuitive": "Statistical Methods in Medical Research", "rue2009approximate": "Journal of the Royal Statistical Society: Series B (Statistical Methodology)",
         "knorrheld2000bayesian": "Statistics in Medicine", "steegen2016increasing": "Perspectives on Psychological Science", "lipsitch2010negative": "Epidemiology",
         "greenland1989ecological": "International Journal of Epidemiology"}

def pandoc_cites(t): return re.sub(r"\s*\[(@[^\]]+)\]", lambda m: " {{" + ",".join(k.strip().lstrip("@") for k in m.group(1).split(";")) + "}}", t)
def person(a):
    p = a.strip().split(" "); return (" ".join(p[:-1]) if len(p) > 1 else p[0]), (p[-1] if len(p) > 1 else "")
def surnames(k):
    au = REFS[k]["authors"]; return [au] if isinstance(au, str) else [person(a)[0] for a in au]
def label(k):
    s = surnames(k); y = REFS[k]["year"]
    return (s[0] if len(s) == 1 else f"{s[0]} and {s[1]}" if len(s) == 2 else f"{s[0]} et al."), y
def aylabel(k):
    a, y = label(k); return ((y, a.lower()), f"{a}, {y}")
def harvard(k):
    r = REFS[k]
    if r["type"] == "web": return f"{r['authors']} ({r['year']}) {r['title']}. Available at: {r['url']} (accessed {r['accessed']})."
    au = r["authors"]; A = ", ".join(au[:3]) + " et al." if len(au) > 3 else (au[0] if len(au) == 1 else ", ".join(au[:-1]) + " and " + au[-1])
    t = r["title"].rstrip("."); s = f"{A} ({r['year']}) {t}{'' if t.endswith(('?', '!')) else '.'} <i>{JFULL.get(k, r.get('journal_full') or r['journal'])}</i>"
    if r.get("vol"): s += f" {r['vol']}" + (f"({r['issue']})" if r.get("issue") else "") + (f": {str(r['pages']).replace('-', chr(8211))}" if r.get("pages") else "")
    return s + "." + (f" DOI: {r['doi']}." if r.get("doi") else "")
def hsort(k):
    au = REFS[k]["authors"]; return ([au.lower()] if isinstance(au, str) else [a.lower() for a in au], REFS[k]["year"])

class Doc(JKMS):
    raw = False
    def __init__(self, **k):
        super().__init__(**k)
        st = self.doc.styles["Normal"]; st.font.size = Pt(10)
        for s in self.doc.sections: s.top_margin = s.bottom_margin = s.left_margin = s.right_margin = Cm(3.0)
    def aylabel(self, k): return aylabel(k)
    def para(self, text, sz=10, **k): return super().para(text, sz=sz, **k)
    def _emit(self, p, text, sz=10, bold=False, italic=False):
        if not self.raw: text = re.sub(r"(^|[\s(;])-(?=\d)", "\\1\u2212", text)
        return super()._emit(p, text, sz, bold, italic)
    def save(self, path):   # CDH: 표 9pt 단일 간격
        from docx.enum.text import WD_LINE_SPACING
        for t in self.doc.tables:
            for r in t.rows:
                for c in r.cells:
                    for p in c.paragraphs: p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        for p in self.doc.paragraphs:
            if p.text.startswith(("Values are means", "Odds ratios (95% credible intervals) per standard deviation from", "Odds ratios (95% credible intervals) per standard deviation.", "Robust:", "Age-sex-standardised prevalence across", "WAIC,", "Manufacturing employment (2017")):
                p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
                for r in p.runs: r.font.size = Pt(9)
        return super().save(path)
    def new_page(self, orient="portrait"):
        s = super().new_page(orient); s.top_margin = s.bottom_margin = s.left_margin = s.right_margin = Cm(3.0); return s

H = C.HYP; FORCED = ["urban_population_share", "clinics_per_100k", "sex_ratio", "fiscal_independence", "pop80_share"]
def f3(d): return f"{C.r3(d['OR'])} ({C.r3(d['lo'])}–{C.r3(d['hi'])})"
def table1():
    P = pd.read_csv(os.path.join(AN, "panel_BRSH_std3.csv"), encoding="utf-8-sig")
    META = pd.read_csv(os.path.join(AN, "panel_META_BRSH_std3.csv"), encoding="utf-8-sig")
    FO = pd.read_csv(os.path.join(ROOT, "02_패널", "BRSH_fh_outcomes.csv"), encoding="utf-8-sig")
    P = P.merge(FO[["year", "fbd_region", "SR", "N"]], on=["year", "fbd_region"], validate="one_to_one")
    win = [("2014–2016", P.year <= 2016), ("2017–2019", P.year.between(2017, 2019)), ("2020–2024", P.year >= 2020), ("2014–2024", P.year > 0)]
    hdr = ["Variable"] + [w for w, _ in win]
    def ms(x): return f"{x.mean():.1f} ({x.std():.1f})"
    rows = [["__H__", "Outcome"], ["Toothbrushing after lunch, %"] + [ms(P.SR[m]) for _, m in win], ["Respondents per district-year"] + [f"{P.N[m].median():.0f}" for _, m in win], ["__H__", "Covariates"]]
    rawcol = {r.variable: r.raw_column if "raw_column" in META else None for r in META.itertuples()}
    for v in H + FORCED:
        col = v + "_raw" if v + "_raw" in P else None
        if col is None: rows.append([C.LABEL[v]] + ["–"] * len(win)); continue
        x = np.exp(P[col]) if v == C.DEN and P[col].max() < 10 else P[col]
        rows.append([("Dentists per 100,000" if v == C.DEN else C.LABEL[v])] + [ms(x[m]) for _, m in win])
    return hdr, rows
def table2():
    hdr = ["Covariate (expected direction)", "Main model", "Subsets credible, %", "Type II", "Type IV", "Between districts", "Within districts"]
    rows = []; MV = N["mv"]
    EXP = {C.W: "+", C.MV: "+", C.FAC: "+", C.SM: "−", C.DR: "−", C.ED: "+", C.LIV: "−", C.DEN: "+", C.MF: "−"}
    for g, vs in C.GROUPS:
        rows.append(["__H__", g])
        for v in vs:
            m = N["main"][v]
            if v in H:
                rows.append([f"{C.LABEL[v]} ({EXP[v]})", f3(m), f"{100*MV[v]['cred']:.0f}", f3(N['kh'][f'KH_II|{v}']), f3(N['kh'][f'KH_IV|{v}']), f3(N['bw'][f'{v}|between']), f3(N['bw'][f'{v}|within'])])
            else:
                rows.append([C.LABEL[v], f3(m), "", f3(N['kh'][f'KH_II|{v}']), f3(N['kh'][f'KH_IV|{v}']), "", ""])
    return hdr, rows
T1N = "Values are means (standard deviations) of district-year values unless stated otherwise. Outcome prevalence is age-sex-standardised. Survey covariates are age-sex-standardised prevalences (%). Values are before standardisation; moderate or vigorous activity (2018) and manufacturing employment (2014–2016) include interpolated or carried values."
T2N = (f"Odds ratios (95% credible intervals) per standard deviation from spatio-temporal Fay–Herriot models with the five controls; {len(N['kh'])//14 if False else 9} covariates were entered together. "
       "Random subsets: share of 200 random covariate subsets containing the covariate in which it was credible. Type II and type IV: Knorr-Held space-time interaction models. "
       "Between and within districts: district means and deviations from them entered separately.")
FIGCAP = [
 "Figure 1. District toothbrushing after lunch among adults in Korea, 2014–2024. (A) Mean age-sex-standardised prevalence. (B) Spatial odds ratios with controls only. (C) Spatial odds ratios in the main model. (D) District-specific trend, odds ratio per 10 years. Solid outlines mark districts with a posterior probability above 0.95 of a positive effect and dashed outlines those below 0.05.",
 "Figure 2. Odds ratios per standard deviation (95% credible intervals) for district covariates in the main model, the Knorr-Held type IV model, between-district and within-district components, and the model adjusted for county status and population density. Filled symbols are credible; † marks covariates meeting the robustness criteria.",
 "Figure 3. Toothbrushing after lunch among districts by survey year, 2014–2024. (A) Mean and 10th and 90th percentiles of district prevalence. (B) Difference and ratio between the 90th and 10th percentiles.",
]

def cites(t): return pandoc_cites(t)
def wc(t): return len(re.sub(r"\[@[^\]]+\]|<[^>]+>", "", t).split())
def body_texts(): return C.INTRO + C.METHOD + C.RESULTS + C.DISCUSSION + [C.CONCLUSION]

def results_fill():
    R = list(C.RESULTS); rsr = N.get("rsr", {}); ro = N.get("rural_OR", {})
    lost = [C.LABEL[v].lower() for v in C.ROB if not rsr.get(v, {}).get("cred", True)]
    rs = (f"With restricted spatial regression, {', '.join(lost)} {'was' if len(lost)==1 else 'were'} no longer credible." if lost else "All four remained credible with restricted spatial regression.") if rsr else ""
    gain = [C.LABEL[v].lower() for v in C.HYP if v not in C.ROB and rsr.get(v, {}).get("cred")]
    if gain: rs += f" Under restricted spatial regression, {', '.join(gain)} became credible ({C.o(rsr[[v for v in C.HYP if C.LABEL[v].lower()==gain[0]][0]])})."
    lostu = [C.LABEL[v].lower() for v in C.ROB if not ro.get(v, {}).get("cred", True)]
    us = (f"After adjustment for county status and population density, {', '.join(lostu)} {'was' if len(lostu)==1 else 'were'} no longer credible." if lostu else "All four remained credible after adjustment for county status and population density.") if ro else ""
    R[2] = R[2].replace("RSR_SENT", rs).replace("URBAN_SENT", us).replace("  ", " ")
    return R

def build_ms(anon=True):
    J = Doc(cite="ay")
    J.para(f"<b>{C.TITLE}</b>")
    J.page_break(); J.head("Abstract", 1)
    for h, t in C.ABSTRACT: J.para(f"<b><i>{h}:</i></b> {t}")
    J.para(f"<b>Keywords:</b> {', '.join(C.KEYWORDS)}", before=6)
    J.page_break(); J.counting = True
    J.head("Introduction", 1)
    for p in C.INTRO: J.para(cites(p), indent=True)
    J.head("Method", 1)
    for p in C.METHOD: J.para(cites(p), indent=True)
    J.head("Results", 1)
    for p in results_fill(): J.para(cites(p), indent=True)
    J.head("Discussion", 1)
    for p in C.DISCUSSION: J.para(cites(p), indent=True)
    J.head("Conclusion", 1); J.para(cites(C.CONCLUSION), indent=True)
    J.counting = False
    keys = sorted(set(J.cites), key=hsort)
    J.page_break(); J.head("References", 1); J.raw = True
    for k in keys: J.para(harvard(k))
    J.raw = False
    J.page_break(); J.head("Figure legends", 1)
    for c in FIGCAP: J.para(c)
    J.new_page("landscape"); h, r = table1(); J.table_fit("Table 1. Toothbrushing after lunch and district covariates, Korea, 2014–2024", h, r, T1N, fs_max=9)
    J.new_page("landscape"); h, r = table2(); J.table_fit("Table 2. District covariates and toothbrushing after lunch: main model and robustness checks", h, r, T2N, fs_max=9)
    return J, keys

def build_title():
    J = Doc(cite="ay")
    J.para(f"<b>{C.TITLE}</b>")
    J.para(f"Running title: {C.SHORT}", before=6)
    J.para("Seongdae Kim", before=12)
    J.para(AFFIL)
    J.para(f"Corresponding author: Seongdae Kim, {ADDR}. Telephone: +82-2-320-1114. Email: {EMAIL}. ORCID: {ORCID}", before=6)
    J.para(f"Word count: {WC} words (main text); abstract {AWC} words; 2 tables; 3 figures; {NREF} references.", before=6)
    J.head("Acknowledgements", 1); J.para("Not applicable.")
    J.head("Author contributions", 1); J.para("Seongdae Kim: conceptualisation, data curation, formal analysis, methodology, software, visualisation, writing – original draft, writing – review and editing.")
    J.head("Statements and Declarations", 1)
    J.head("Ethical considerations", 2); J.para("The Institutional Review Board of Hongik University approved the study (approval number 7002340-202511-HR-012; continuing review 7002340-202611-HR-040) and waived informed consent because only published aggregate statistics were used.")
    J.head("Consent to participate", 2); J.para("Not applicable. The study used only published aggregate statistics; the requirement for informed consent was waived by the Institutional Review Board of Hongik University.")
    J.head("Consent for publication", 2); J.para("Not applicable.")
    J.head("Declaration of conflicting interest", 2); J.para("The author declared no potential conflicts of interest with respect to the research, authorship, and/or publication of this article.")
    J.head("Funding statement", 2); J.para("The author received no financial support for the research, authorship, and/or publication of this article.")
    J.head("Data availability", 2); J.para(f"The data are publicly available from the Korean Statistical Information Service (https://kosis.kr). The analysis code is available at {GITHUB} and archived at https://doi.org/{ZENODO}. No individual-level data were used.")
    J.head("Use of artificial intelligence", 2); J.para("Claude (Anthropic) and ChatGPT (OpenAI) were used only for translation and English-language editing; the author reviewed all output and takes full responsibility for the content.")
    return J

def pdf(path):
    subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", os.path.dirname(path), path], capture_output=True, timeout=300)

if __name__ == "__main__":
    import sanitize_docx
    J, keys = build_ms()
    WC = sum(wc(t) for t in body_texts()); AWC = sum(len(f"{h}: {t}".split()) for h, t in C.ABSTRACT); NREF = len(keys)
    fm = os.path.join(OUT, "Manuscript.docx"); J.save(fm)
    T = build_title(); ft = os.path.join(OUT, "Title_page.docx"); T.save(ft)
    open(os.path.join(OUT, "wordcount.txt"), "w").write(f"main {WC}\nabstract {AWC}\nrefs {NREF}\n")
    open(os.path.join(OUT, "reference_list.txt"), "w").write("\n".join(re.sub(r"<[^>]+>", "", harvard(k)) for k in keys))
    for f in (fm, ft): pdf(f)
    sanitize_docx.sanitize(fm, anon=True); sanitize_docx.sanitize(ft, anon=False)
    print(f"[build_cdh] main {WC} words, abstract {AWC}, refs {NREF}: {keys}")
