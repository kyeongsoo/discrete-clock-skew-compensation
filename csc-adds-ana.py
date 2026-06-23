#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file     csc-adds-ana.py
# @author   Kyeong Soo (Joseph) Kim <kyeongsoo.kim@gmail.com>
# @date     2026-04-23
#
# @brief    Plot the results of CSC ADDS analysis.
#
# @remarks The results are published in the following paper:
# - Kyeong Soo Kim, "Space-time trade-off in integer linear scaling
#   rounded to the nearest integer through multiplicative and additive
#   decomposition," arXiv e-prints arXiv:2605.21400v [cs.DS], May 2026.
#   [Online]. Available: https://arxiv.org/abs/2605.21400
#
# @remarks  Copyright (C) 2026 Kyeong Soo (Joseph) Kim. All rights reserved.
#
# @remarks  SPDX-License-Identifier: MIT
#

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc

# global constant
USE_TEX = True  # set to 'False' if LaTeX is not available

# ----------------------------------
# post-process and plot the results
# ----------------------------------
if USE_TEX:
    rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    rc("text", usetex=True)
rc("lines", linewidth=0.5, markersize=3)

# dfs = []
xticklabels = [r"$10^{0}$", r"$10^{1}$", r"$10^{2}$", r"$10^{3}$"]
for int_type in ["int32", "int64"]:
    files = list(glob.glob(f"./out/csc-adds-ana_{int_type}_*.bin"))
    csc_errs = np.zeros((1000, len(files)), dtype=getattr(np, int_type))
    for i, file in enumerate(files):
        csc_errs[:, i] = np.fromfile(file, dtype=getattr(np, int_type))
    # df = pd.DataFrame(csc_errs, columns=col_names[0] if int_type == 'int32' else col_names[1])
    df = pd.DataFrame(csc_errs)

    ax = df.boxplot(showfliers=False)
    ax.set_xticklabels(xticklabels)
    # ax.set_yscale('log')
    plt.xlabel(r"$N$")
    plt.ylabel("Integer Linear Scaling Error")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.savefig(f"./out/csc-adds-ana_{int_type}.pdf")
    plt.clf()
