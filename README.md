Minimal Statistical Reproducibility
Status: LOCAL AUTHOR-REVIEW CANDIDATE. NOT APPROVED FOR PUBLIC RELEASE.

Scope:
- Recalculate patient-level Table 2 method medians and paired effects.
- Recalculate 20 E7.1 Wilcoxon/Holm comparisons.
- Recalculate six Table 2 percentile bootstrap 95% confidence intervals.
- Check selected E6 calibration perturbation and repetition results.
- Recalculate eight E6 matched temporal effects, bootstrap intervals,
  and Holm-adjusted p-values.

Run:
    python verify_statistics.py

Tested analysis environment:
- Python 3.10.20
- NumPy 2.2.6
- SciPy 1.15.3

Data:
- Input tables are frozen derived project results, not raw gait data.
- Patient identifiers in this candidate are consistently replaced
  with local P001–P050 aliases.
- These aliases are not a guarantee of complete anonymization.

Important limitations:
- Does not recompute marker reconstruction or C3 registration from C3D.
- Does not independently reproduce the original E2 or E3 pipelines.
- Does not include C3 models, GLB assets, Qt software, or raw StrokePiG data.
- Does not reproduce all E6 perturbation confidence intervals from
  the underlying random pose perturbations.
- Reproduces selected statistical summaries from supplied derived data.

Original dataset:
Van Criekinge et al., Scientific Data, 2023.
https://doi.org/10.1038/s41597-023-02767-y
Stroke C3D record:
https://doi.org/10.6084/m9.figshare.24192483.v1

The candidate must be reviewed for sharing permission,
patient-level disclosure risk, source provenance and completeness
before any publication.
