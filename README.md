# Hybrid RPRG-Initialized NCME for Nonlinear Chirp Mode Estimation

This repository contains a MATLAB implementation of a hybrid pipeline for estimating and reconstructing multi-component **nonlinear chirp** modes under **crossings** and **large dynamic range**, where weak modes are often buried by residual artifacts from strong modes.

**Core idea:** use **RPRG** for crossing-aware ridge initialization, then apply **NCME (ADMM-style multi-mode joint refinement)** to jointly refine IF trajectories and reconstruct modes, improving weak-mode recovery.

---

## Method Overview

Main function is hybrid_RPRG_VNCMD_NCME

Pipeline (high level):

1. Compute time–frequency representation (STFT)
2. Extract strong-mode ridges via sequential ridge extraction, then **RPRG** to handle crossings
3. (Optional) refine ridge initialization with VNCMD
4. Use **NCME_multi** on the original signal to jointly refine strong modes
5. Add weak modes incrementally: detect ridge on the residual, then re-run **NCME_multi** for joint refinement

---

## Results (Examples)

### (A) Strong/weak separation + reconstruction
run("test6.m")
![Strong/weak separation](figures/fig_strongweak_combo.png)

### (B) Crossing case: regrouping + joint refinement
run("Test3.m")

![Crossing](figures/fig_cross.png)

### (C) Multiple crossings: refined IF tracking
run("testncme2.m")
![Multi-crossing](figures/fig_multicross.png)



## Quick Start

### Requirements
- MATLAB R20XXa (tested on: MATLAB R2025b (Windows 11))
- Dependencies: `STFT`, `extridge_mult`, `RPRG`, `NCME_multi`, `findridges`, `Dechirp_filter`, `curvesmooth`, (optional) `VNCMD`

> **Note:** Some dependencies are external implementations. See **Acknowledgements** below.


