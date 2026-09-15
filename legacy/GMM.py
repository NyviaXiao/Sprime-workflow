import os
import glob
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import chi2
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize

def format_pvalue(p):
    if p < 1e-4:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"

# =========================
# 1. 拟合 GMM 并返回 log-likelihood 和参数
# =========================
def fit_gmm(xdata, n_components, random_state=0):
    X = xdata.reshape(-1, 1)
    n = X.shape[0]

    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        random_state=random_state
    )
    gmm.fit(X)

    logL = gmm.score(X) * n

    params = {
        "weights": gmm.weights_,
        "means": gmm.means_.flatten(),
        "stds": np.sqrt(gmm.covariances_.flatten())
    }

    return logL, params


# =========================
# 2. 对 1 / 2 / 3 pulse 做逐级 LRT
# =========================
def fit_gmm_lrt_123(xdata):
    # 拟合 1, 2, 3 个成分
    logL1, p1 = fit_gmm(xdata, 1)
    logL2, p2 = fit_gmm(xdata, 2)
    logL3, p3 = fit_gmm(xdata, 3)

    # ---------- 1 vs 2 ----------
    LR12 = 2 * (logL2 - logL1)
    df12 = 3
    p_value_12 = 1 - chi2.cdf(LR12, df12)

    # Bonferroni（和你之前一样）
    alpha12 = 0.05 / 23

    if p_value_12 < alpha12:
        # ---------- 2 vs 3 ----------
        LR23 = 2 * (logL3 - logL2)
        df23 = 3
        p_value_23 = 1 - chi2.cdf(LR23, df23)
        alpha23 = 0.05 / 46

        if p_value_23 < alpha23:
            n_components = 3
            params = p3
        else:
            n_components = 2
            params = p2
    else:
        # ---------- 1 vs 3 ----------
        LR13 = 2 * (logL3 - logL1)
        df13 = 6
        p_value_23 = 1 - chi2.cdf(LR13, df13)
        alpha23 = 0.05 / 46

        if p_value_23 < alpha23:
            n_components = 3
            params = p3
        else:
            n_components = 1
            params = p1

    return {
        "n_components": n_components,
        "weights": params["weights"],
        "means": params["means"],
        "stds": params["stds"],
        "p_value_12": p_value_12,
        "p_value_23": p_value_23
    }


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# =========================
# 目标群体
# =========================
target_groups = ["CMBR", "LALL", "LALS", "THW2"]

# =========================
# 创建 2×2 图
# =========================
fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=600)
axes = axes.flatten()

# =========================
# 统一 x 轴范围（推荐）
# =========================
xmin, xmax = 0, 1.0

# =========================
# 主循环
# =========================
for ax, group in zip(axes, target_groups):

    filepath = f"{group}_19_xiao_twodeni.txt"
    data = pd.read_table(filepath, sep="\t")

    # ---------- 数据过滤 ----------
    xdata = data[
        ((data.iloc[:, 5] > 0.3) | (data.iloc[:, 6] > 0.3)) &
        (data.iloc[:, 4] < 0.3)
    ].iloc[:, 5].to_numpy()

    if len(xdata) < 10:
        ax.set_visible(False)
        continue

    # ---------- LRT ----------
    lrt_result = fit_gmm_lrt_123(xdata)
    k = lrt_result["n_components"]
    p12 = lrt_result["p_value_12"]
    p23 = lrt_result["p_value_23"]

    # ---------- GMM ----------
    X = xdata.reshape(-1, 1)
    gmm = GaussianMixture(n_components=k, random_state=0)
    gmm.fit(X)

    weights = gmm.weights_
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())

    # ---------- 直方图 ----------
    ax.hist(
        xdata,
        bins=25,
        density=True,
        alpha=0.5,
        edgecolor="black",
        linewidth=0.5
    )

    # ---------- GMM 曲线 ----------
    x_grid = np.linspace(xmin, xmax, 500)
    mixture_pdf = np.zeros_like(x_grid)

    colors = ["#4C72B0", "#DD8452", "#55A868"]  # 分量颜色

    for i, (w, mu, sd) in enumerate(zip(weights, means, stds)):
        comp_pdf = w * norm.pdf(x_grid, mu, sd)
        mixture_pdf += comp_pdf
        ax.plot(x_grid, comp_pdf, linestyle="--", linewidth=1.5, color=colors[i])

    # 主分布
    ax.plot(x_grid, mixture_pdf, linewidth=2.5, color="black")

    # ---------- 标题 ----------
    ax.set_title(f"{group}  (k = {k})", fontsize=14)

    # ---------- p-value ----------
    text_lines = [f"1 vs 2: {format_pvalue(p12)}"]
    if k >= 2:
        text_lines.append(f"2 vs 3: {format_pvalue(p23)}")

    ax.text(
        0.65, 0.92,
        "\n".join(text_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11
    )

    # ---------- 坐标轴 ----------
    ax.set_xlim(xmin, xmax)

    # 四边框
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

# =========================
# 全局标签
# =========================
fig.text(0.5, 0.04, "Match rate to Denisovan", ha="center", fontsize=14)
fig.text(0.04, 0.5, "Density", va="center", rotation="vertical", fontsize=14)

# =========================
# 布局
# =========================
plt.tight_layout(rect=[0.05, 0.05, 1, 1])
plt.show()