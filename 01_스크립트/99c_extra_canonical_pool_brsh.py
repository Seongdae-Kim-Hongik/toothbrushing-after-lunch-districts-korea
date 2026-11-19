#!/usr/bin/env python3
#  tertiary_educated_share = 대졸 이상 비율(같은 파일)
#                            2014-16 원표 없음 → 각 시군구 2017 값으로 고정(뒤로 이월, 규칙 기록), 결측 시군구는 같은 해 전국 중앙값
#  brsh_SR                 = 결과(02_패널/BRSH_fh_outcomes.csv SR, 세종 원표 복원) — 80의 outcome_column 용
#  산출: 02_패널/EXT_canonical_pool_BRSH_261119.csv (2011-2024 격자; 결과는 2014+만 값)
import os, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..")); PAN = os.path.join(ROOT, "02_패널")
rd = lambda f: (lambda x: x.assign(fbd_region=x.fbd_region.map(NFC)))(pd.read_csv(f, encoding="utf-8-sig"))
O = rd(os.path.join(PAN, "EXT_canonical_pool_OBES_261111.csv")); D = rd(os.path.join(PAN, "EXT_canonical_pool_PERIO_261109.csv"))
MF = rd(os.path.join(ROOT, "04_분석", "SCRN", "EXT_manufacturing_share_2017_24.csv")); B = rd(os.path.join(PAN, "BRSH_fh_outcomes.csv"))
out = O[["year", "fbd_region", "basic_livelihood_rate", "tertiary_educated_share"]].merge(D[["year", "fbd_region", "dentists_per100k"]], on=["year", "fbd_region"], how="inner", validate="one_to_one")
assert len(out) == 229 * 14 and out.fbd_region.nunique() == 229
out = out.merge(MF, on=["year", "fbd_region"], how="left", validate="one_to_one")
n0 = int(out.loc[out.year >= 2017, "manufacturing_share"].isna().sum()); print(f"[99c-brsh] 제조업 비율 2017-24 결측 {n0}칸:", out.loc[(out.year >= 2017) & out.manufacturing_share.isna(), ["year", "fbd_region"]].values.tolist())
out["manufacturing_share"] = out.manufacturing_share.fillna(out.groupby("year").manufacturing_share.transform("median")).where(out.year >= 2017)
v17 = out[out.year == 2017].set_index("fbd_region").manufacturing_share
out.loc[out.year < 2017, "manufacturing_share"] = out.loc[out.year < 2017, "fbd_region"].map(v17)
out = out.merge(B[["year", "fbd_region", "SR"]].rename(columns={"SR": "brsh_SR"}), on=["year", "fbd_region"], how="left", validate="one_to_one")
assert out.loc[out.year >= 2014, "brsh_SR"].notna().all() and out.loc[out.year < 2014, "brsh_SR"].isna().all()
ref = out[out.year.isin([2020, 2024])].groupby("fbd_region").basic_livelihood_rate.mean(); yy = out.year.isin([2021, 2022, 2023])
r = out.loc[yy, "basic_livelihood_rate"] / out.loc[yy, "fbd_region"].map(ref); assert r.between(0.2, 5).all(); print(f"[99c-brsh] 수급률 2021-23 비 범위 {r.min():.2f}-{r.max():.2f}")
for v in ["dentists_per100k", "manufacturing_share"]:
    print(f"[99c-brsh] 게이트 {v}: 0이하 {int((out[v] <= 0).sum())}, min {out[v].min():.2f} max {out[v].max():.2f}")
cols = ["dentists_per100k", "basic_livelihood_rate", "tertiary_educated_share", "manufacturing_share"]
assert out[cols].notna().all().all()
out = out[["year", "fbd_region"] + cols + ["brsh_SR"]].sort_values(["fbd_region", "year"])
fp = os.path.join(PAN, "EXT_canonical_pool_BRSH_261119.csv"); out.to_csv(fp, index=False, encoding="utf-8-sig"); print("[99c-brsh] 저장", fp, out.shape)
print(out[out.year >= 2014][cols].describe().T[["min", "50%", "max"]].round(3))
