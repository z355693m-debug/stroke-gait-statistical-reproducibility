# Stroke Gait Statistical Reproducibility

This repository contains the minimal statistical reproducibility materials for the study **“Patient-specific geometry in stroke gait: Cross-trial validation and calibration sensitivity.”**

The study retrospectively analyzes gait recordings from 50 people with stroke, covering 149 dynamic gait trials and 298 within-patient directed calibration-to-held-out pairs. The files here are intended to verify **selected patient-level statistics from frozen derived results**, not to reproduce the complete motion-capture-to-model pipeline.

## Contents

```text
stroke-gait-statistical-reproducibility/
├── README.md
├── LICENSE
├── verify_statistics.py
└── data/
    ├── LICENSE.md
    ├── e1.csv
    ├── e4.csv
    ├── information_paired.csv
    ├── information_repeats.csv
    ├── multiplicity.csv
    ├── orientation.csv
    ├── orientation_repeats.csv
    └── position.csv
```

The eight CSV files are derived statistical tables, not raw gait recordings. Patient identifiers have been replaced by consistent local aliases (`P001`–`P050`). Pseudonymization **does not guarantee complete anonymization**; reuse must respect applicable data protection requirements.

## What can be verified

Run `verify_statistics.py` to:

1. Recalculate the Table 2 patient-level method medians and paired improvements for the eight-marker reconstruction (E1) and model-defined C3 anchor localization (E4) endpoints.
2. Recalculate 20 E7.1 Wilcoxon comparisons and their Holm-adjusted p-values.
3. Recalculate six patient-bootstrap 95% confidence intervals reported in Table 2.
4. Check selected E6 position/orientation perturbation results and repetition-level heterogeneity.
5. Recalculate eight E6 matched temporal-subsampling effects, patient-bootstrap confidence intervals, and Holm-adjusted p-values.

Statistical inference is performed at the **patient level**. Trials, directed pairs, frames, channels, and perturbation repetitions are not treated as additional independent patients.

## Requirements

The verification workflow has been tested with:

- **Ubuntu:** Python 3.10.20
- **Windows:** Python 3.11.9
- **NumPy:** 2.2.6
- **SciPy:** 1.15.3

No original laboratory software or external model assets are needed for the statistical verification provided here.

## Run the verification

Download the repository (or clone it), open a terminal in its root directory, and install the specified dependencies:

```bash
python -m pip install numpy==2.2.6 scipy==1.15.3
python verify_statistics.py
```

Using an isolated Python virtual environment is recommended. On systems where `python` points to a different installation, use the corresponding Python executable.

A successful run ends with:

```text
TABLE 2 METHODS/PAIRED EFFECTS: PASS
E7 WILCOXON/HOLM: 20/20 PASS
TABLE 2 BOOTSTRAP 95%CI: 6/6 PASS
E6 PERTURBATION/REPETITION CHECK: PASS
E6 TEMPORAL EFFECT/CI/HOLM: 8/8 PASS
ALL STATISTICAL CHECKS PASSED
```

The verification checks predefined numerical results against the supplied derived CSV files. A `PASS` does not independently establish the correctness of the original experimental data or the complete upstream analysis.

## Source data and attribution

The study uses the publicly available **StrokePiG** gait dataset described by:

Van Criekinge T, et al. (2023). *A full-body motion capture gait dataset of 138 able-bodied adults across the life span and 50 stroke survivors.* **Scientific Data**, 10, 852. https://doi.org/10.1038/s41597-023-02767-y

- Original motion-capture data: https://doi.org/10.6084/m9.figshare.24192483.v1
- Source dataset collection: https://doi.org/10.6084/m9.figshare.c.6503791.v1

The original StrokePiG recordings are **not redistributed** in this repository. When reusing the derived tables, please acknowledge this study and cite the original StrokePiG dataset publication, as appropriate.

## Scope and limitations

This is a **statistical reproducibility package based on supplied derived results**, not a full end-to-end reproduction of the research workflow. In particular, it:

- Does **not** recalculate eight-marker reconstruction or C3 anchor localization from raw C3D recordings.
- Does **not** independently reproduce the upstream E2 offline kinematic-readback or E3 pelvic-width analyses.
- Does **not** include raw StrokePiG records, C3/GLB models, Qt/WebGL applications, or the complete geometry-estimation codebase.
- Does **not** regenerate all E6 perturbation outcomes from the underlying randomly perturbed poses.
- Does **not** demonstrate physical exoskeleton fitting accuracy, anatomical localization accuracy, or actual browser-runtime measurements.

The provided script and CSV files support independently running the **specified statistical checks** on derived data.

## Licensing

- **Statistical verification code:** MIT License; see [`LICENSE`](LICENSE).
- **Derived CSV tables in `data/`:** Creative Commons Attribution 4.0 International (CC BY 4.0); see [`data/LICENSE.md`](data/LICENSE.md).
- **Original StrokePiG dataset:** Governed by the terms of its original repository record; those terms are not changed by this repository.

These licenses apply only to materials for which the contributors hold the necessary rights. Third-party raw datasets and software assets are not licensed or redistributed here.

## Citation

If you use this repository, please cite the associated study once its final bibliographic record is available, and cite the original StrokePiG dataset publication above where relevant. No journal acceptance, publication DOI, or archived software DOI is claimed here.
