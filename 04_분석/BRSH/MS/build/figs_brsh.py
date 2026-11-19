#  그림 규칙: 컬러, 패널 문자 대문자·같은 크기, 600 dpi TIFF, 그림 안 설명문 금지
import os, shutil
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import content_brsh as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(OUT, exist_ok=True)
AN = os.path.normpath(os.path.join(HERE, "..", "..")); N = C.N
im = Image.open(os.path.join(AN, "figs", "Fig_maps_BRSH.tiff")); im.save(os.path.join(OUT, "Figure1.tiff"), compression="tiff_lzw", dpi=(600, 600))
shutil.copy(os.path.join(AN, "figs", "Fig_maps_BRSH.png"), os.path.join(OUT, "Figure1.png"))
plt.rcParams.update({"font.family": "Arial", "font.size": 11})
def save(fig, k):
    for ext, dpi in (("png", 300), ("tiff", 600)): fig.savefig(os.path.join(OUT, f"Figure{k}.{ext}"), dpi=dpi, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None, bbox_inches="tight")
# Figure 2
VARS = C.HYP
SER = [("Main model", lambda v: N["main"][v], "#1f4e79", "o"), ("Knorr-Held type IV", lambda v: N["kh"][f"KH_IV|{v}"], "#2e8b57", "s"),
       ("Between districts", lambda v: N["bw"][f"{v}|between"], "#b8860b", "^"), ("Within districts", lambda v: N["bw"][f"{v}|within"], "#8b1a1a", "v"),
       ("Adjusted for urbanisation", lambda v: N["rural_OR"][v], "#6a3d9a", "D")]
fig, ax = plt.subplots(figsize=(8.5, 9.5)); y0 = np.arange(len(VARS))[::-1]; off = np.linspace(-0.32, 0.32, len(SER))
for (lab, f, col, mk), dy in zip(SER, off):
    for v, y in zip(VARS, y0):
        d = f(v); ax.plot([d["lo"], d["hi"]], [y + dy, y + dy], color=col, lw=1.4)
        ax.plot(d["OR"], y + dy, marker=mk, color=col, mfc=col if d["cred"] else "white", ms=6.5, label=lab if v == VARS[0] else None)
ax.axvline(1, color="#555", lw=0.8, ls="--")
ax.set_yticks(y0); ax.set_yticklabels([C.LABEL[v] + (" †" if v in N["robust"] else "") for v in VARS]); ax.set_xscale("log")
ticks = [0.90, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20, 1.25]; ax.set_xticks(ticks); ax.set_xticklabels([f"{t:.2f}" for t in ticks]); ax.minorticks_off()
ax.set_xlabel("Odds ratio per standard deviation (log scale)"); ax.set_ylim(-0.7, len(VARS) - 0.3)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.42, -0.2), ncol=2, frameon=False, fontsize=10)
plt.tight_layout(); save(fig, 2); plt.close(fig)
# Figure 3
G = pd.DataFrame(N["gap"]).T; G.index = G.index.astype(int); mean = pd.Series({int(k): v for k, v in N["sr"].items()})
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8))
a.fill_between(G.index, G.P10, G.P90, color="#f4a582", alpha=0.45, lw=0, label="10th-90th percentile")
a.plot(G.index, G.P90, color="#b2182b", lw=1.6, marker="o", ms=4, label="90th percentile"); a.plot(mean.index, mean, color="#333", lw=2, marker="o", ms=4, label="Mean")
a.plot(G.index, G.P10, color="#2166ac", lw=1.6, marker="o", ms=4, label="10th percentile")
a.set_ylabel("Toothbrushing after lunch, %"); a.set_ylim(30, 90); a.axvline(2019.5, color="#777", lw=0.8, ls=":"); a.legend(frameon=False, fontsize=10, loc="lower right")
b.plot(G.index, G["diff"], color="#b2182b", lw=2, marker="o", ms=4, label="Difference, percentage points"); b.set_ylabel("90th minus 10th percentile, percentage points"); b.set_ylim(0, 25); b.axvline(2019.5, color="#777", lw=0.8, ls=":")
b2 = b.twinx(); b2.plot(G.index, G.ratio, color="#2166ac", lw=2, marker="s", ms=4, label="Ratio"); b2.set_ylabel("90th to 10th percentile ratio"); b2.set_ylim(1.0, 1.6)
h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels(); b.legend(h1 + h2, l1 + l2, frameon=False, fontsize=10, loc="lower left")
for x in (a, b, b2):
    x.set_xticks(range(2014, 2025, 2)); x.spines["top"].set_visible(False)
for x in (a, b): x.set_xlabel("Survey year")
plt.tight_layout(w_pad=4); fig.canvas.draw()
for lab, x in (("A", a), ("B", b)):
    bb = x.get_position(); fig.text(bb.x0 - 0.06, bb.y1 + 0.02, lab, fontsize=24, fontweight="bold", ha="left", va="bottom")
save(fig, 3); plt.close(fig)
print("figs ok")
