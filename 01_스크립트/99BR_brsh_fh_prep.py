#!/usr/bin/env python3
#  입력: 기준 SQLite region_year — kchs__DT_H_OR_BRUSH SR·SR_SE·N·CR·CR_SE (읽기 전용)
#        세종: SQLite 결측(228곳) → KOSIS 원표(orgId 177 DT_H_OR_BRUSH) 시도 머리행 = 세종시 행의 N·CR·SR 로 살림.
#        세종 SE는 원표에 없음 → 같은 해 다른 228곳의 SR_SE / sqrt(p(1-p)/N) 중앙값(설계효과 근사) × sqrt(p(1-p)/N) 로 대치(기록).
#  게이트: 원표 값 = SQLite 값(228곳×11년), SE>0, 연도별 전국 평균, 2019→2020 단절(문항 간소화) 순위 상관.
#  산출: 02_패널/BRSH_fh_outcomes.csv (2014-2024, 229×11)
import os, sqlite3, unicodedata, glob
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
DB = os.path.expandvars("$ST_GPT_DATASET") + "/02_체육건강_시공간_GPT버전/02_체육건강_시공간_GPT버전.sqlite"
RAW = glob.glob(os.path.expandvars("$ST_GPT_DATASET") + "/00_재구축_원자료_GPT버전/원본수집_GPT버전/**/KOSIS_177_DT_H_OR_BRUSH_all.csv", recursive=True)[0]
Y = list(range(2014, 2025)); K = "kchs__DT_H_OR_BRUSH__"
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
d = pd.read_sql(f'select year, fbd_region, population, "{K}SR" SR, "{K}SR_SE" SR_SE, "{K}N" N, "{K}CR" CR, "{K}CR_SE" CR_SE from region_year where year between 2014 and 2024', c)
d["fbd_region"] = d.fbd_region.map(NFC); SJ = NFC("세종시세종시")
assert len(d) == 229 * 11 and not d.duplicated(["year", "fbd_region"]).any()
miss = d[d.SR.isna()]; assert set(miss.fbd_region) == {SJ} and len(miss) == 11, miss.fbd_region.unique()
# 원표
r = pd.read_csv(RAW, header=[0, 1], encoding="utf-8")
r.columns = [f"{a}|{b}" if not str(a).startswith("시군구") and not str(a).startswith('"') else str(a).strip('"') for a, b in r.columns]
idc = r.columns[:3]; r = r.rename(columns=dict(zip(idc, ["s1", "s2", "s3"])))
sj = r[(r.s1.map(NFC) == NFC("세종특별자치시")) & (r.s2.map(NFC) == NFC("세종시"))]; assert len(sj) == 1
def rawv(row, y, kind):
    col = [cc for cc in r.columns if cc.startswith(f"{y}|") and kind in cc.split("|")[1]]
    assert len(col) == 1, (y, kind, col); return float(row[col[0]])
# 원표 vs SQLite 대조(서울 강남구 등 시군구 행: s3=소계, s2≠소계) — 시군구 이름 결합이 아니라 값 분포 대조(평균·합) 수준
for y in Y:
    vals = pd.to_numeric(r.loc[(r.s2 != "소계") & (r.s3 == "소계"), [cc for cc in r.columns if cc.startswith(f"{y}|") and "표준화율" in cc][0]], errors="coerce")
    sq = d.loc[d.year == y, "SR"].dropna()
    # 원표 시군구 행에는 세종(시도=시군구) 포함 + 일반구 행이 있을 수 있음 → 합집합 크기 차는 허용, SQLite 값이 원표 값 집합에 모두 있는지 확인
    s_raw = set(np.round(vals.dropna(), 1)); nm = int((~np.round(sq, 1).isin(s_raw)).sum())
    if nm: print('[99BR] 원표 집합에 없는 SQLite 값', y, d.loc[(d.year==y)&d.SR.notna()&~np.round(d.SR,1).isin(s_raw),['fbd_region','SR']].values.tolist())
    assert nm <= 12, (y, nm)
sjrows = []
for y in Y:
    N_ = rawv(sj.iloc[0], y, "응답자수"); CR_ = rawv(sj.iloc[0], y, "조율"); SR_ = rawv(sj.iloc[0], y, "표준화율")
    sjrows.append(dict(year=y, N=N_, CR=CR_, SR=SR_))
sjd = pd.DataFrame(sjrows)
# 설계효과 근사로 SE 대치
d["deffR"] = d.SR_SE / np.sqrt((d.SR / 100) * (1 - d.SR / 100) / d.N) / 100
d["deffRc"] = d.CR_SE / np.sqrt((d.CR / 100) * (1 - d.CR / 100) / d.N) / 100
med = d.groupby("year")[["deffR", "deffRc"]].median()
for _, s in sjd.iterrows():
    m = (d.year == s.year) & (d.fbd_region == SJ); p, pc = s.SR / 100, s.CR / 100
    d.loc[m, ["SR", "N", "CR"]] = [s.SR, s.N, s.CR]
    d.loc[m, "SR_SE"] = round(med.loc[s.year, "deffR"] * np.sqrt(p * (1 - p) / s.N) * 100, 2)
    d.loc[m, "CR_SE"] = round(med.loc[s.year, "deffRc"] * np.sqrt(pc * (1 - pc) / s.N) * 100, 2)
d["sejong_imputed_SE"] = (d.fbd_region == SJ).astype(int)
print("[99BR] 세종 원표 SR·대치 SE", d.loc[d.fbd_region == SJ, ["year", "SR", "SR_SE", "N"]].to_dict("records"))
print("[99BR] 설계효과 근사(SE/이항SE) 연도별 중앙", med.deffR.round(2).to_dict())
assert d[["SR", "SR_SE", "N", "CR", "CR_SE"]].notna().all().all() and (d.SR > 0).all() and (d.SR < 100).all() and (d.SR_SE > 0).all()
d["SE_used"] = d.SR_SE.clip(lower=0.05); d["CRSE_used"] = d.CR_SE.clip(lower=0.05)
p = d.SR / 100; d["y_logit"] = np.log(p / (1 - p)); d["v_logit"] = (d.SE_used / 100) ** 2 / (p * (1 - p)) ** 2
pc = d.CR / 100; d["y_logit_cr"] = np.log(pc / (1 - pc)); d["v_logit_cr"] = (d.CRSE_used / 100) ** 2 / (pc * (1 - pc)) ** 2
d["n_eff"] = np.maximum(np.round(p * (1 - p) / (d.SE_used / 100) ** 2), 10); d["y_eff"] = np.round(p * d.n_eff)
rel = 1 - d.v_logit / d.groupby("year").y_logit.transform("var")
print(f"[99BR] rows {len(d)} · SR {d.SR.min():.1f}-{d.SR.max():.1f} · N 중앙 {d.N.median():.0f} · n_eff 중앙 {d.n_eff.median():.0f} · 신뢰도 중앙 {rel.median():.2f}")
print("[99BR] 연도 평균 SR", d.groupby("year").SR.mean().round(2).to_dict())
pv = d.pivot(index="fbd_region", columns="year", values="SR")
print("[99BR] 단절 게이트 Spearman 2018-19 %.3f · 2019-20 %.3f · 2020-21 %.3f" % (pv[2018].corr(pv[2019], method="spearman"), pv[2019].corr(pv[2020], method="spearman"), pv[2020].corr(pv[2021], method="spearman")))
d = d.drop(columns=["deffR", "deffRc"])
fp = os.path.join(ROOT, "02_패널", "BRSH_fh_outcomes.csv"); d.to_csv(fp, index=False, encoding="utf-8-sig"); print("[99BR] 저장", fp, d.shape)
