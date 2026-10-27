
from pathlib import Path
from collections import defaultdict
import csv
import numpy as np
import scipy
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

def load(name):
    with (DATA / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def median(values):
    return float(np.median(np.asarray(values, dtype=float)))

def close(a, b, tol=0.005):
    return abs(float(a) - float(b)) <= tol

def pclose(a, b):
    return abs(float(a)-float(b)) <= max(
        1e-18, abs(float(b))*1e-5
    )

def require(condition, label):
    if not condition:
        raise RuntimeError("FAIL: " + label)

def signed_p(values):
    return float(wilcoxon(
        np.asarray(values, dtype=float),
        alternative="two-sided",
        zero_method="wilcox"
    ).pvalue)

def holm(values):
    p = np.asarray(values, dtype=float)
    result = np.zeros(len(p))
    previous = 0.0
    for rank, index in enumerate(
        np.argsort(p, kind="stable")
    ):
        previous = max(
            previous, min(1.0, (len(p)-rank)*p[index])
        )
        result[index] = previous
    return result

def bootstrap(values, seed):
    a = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    index = rng.integers(
        0, len(a), size=(20000, len(a))
    )
    return np.percentile(
        np.median(a[index], axis=1),
        [2.5, 97.5]
    )

print("PYTHON NUMPY SCIPY:", np.__version__, scipy.__version__)

TABLE = {
    "E1": {
        "methods": (17.02, 30.09, 50.91, 38.61),
        "effects": (9.70, 29.65, 23.03),
        "positive": (47, 50, 50),
    },
    "E4": {
        "methods": (149.51, 154.11, 157.95, 156.95),
        "effects": (5.16, 6.43, 7.24),
        "positive": (46, 49, 48),
    }
}

METHODS = (
    "GAIT_CALIBRATED",
    "STATIC_ONLY",
    "LOPO_GENERIC",
    "LOPO_ANTHRO",
)

COLUMNS = {
    ("GAIT_CALIBRATED", "STATIC_ONLY"):
        "GAIT_vs_STATIC_improvement_mm",
    ("GAIT_CALIBRATED", "LOPO_GENERIC"):
        "GAIT_vs_LOPO_GENERIC_improvement_mm",
    ("GAIT_CALIBRATED", "LOPO_ANTHRO"):
        "GAIT_vs_LOPO_ANTHRO_improvement_mm",
    ("STATIC_ONLY", "LOPO_GENERIC"):
        "STATIC_vs_LOPO_GENERIC_improvement_mm",
    ("STATIC_ONLY", "LOPO_ANTHRO"):
        "STATIC_vs_LOPO_ANTHRO_improvement_mm",
    ("LOPO_GENERIC", "LOPO_ANTHRO"):
        "LOPO_GENERIC_vs_LOPO_ANTHRO_improvement_mm",
}

patient_rows = {
    "E1": load("e1.csv"),
    "E4": load("e4.csv"),
}

def effects(stage, better, comparator):
    pair = (better, comparator)
    reverse = (comparator, better)
    if pair in COLUMNS:
        col, sign = COLUMNS[pair], 1
    elif reverse in COLUMNS:
        col, sign = COLUMNS[reverse], -1
    else:
        raise RuntimeError("Unknown contrast " + str(pair))
    return np.asarray([
        sign * float(row[col])
        for row in patient_rows[stage]
    ])

# A. Table 2 medians, paired effects and positive patients.
for stage, rows in patient_rows.items():
    require(len(rows) == 50, stage + " patient count")
    require(
        len({r["patient"] for r in rows}) == 50,
        stage + " unique patients"
    )
    for i, method in enumerate(METHODS):
        require(
            close(
                median([r[method] for r in rows]),
                TABLE[stage]["methods"][i]
            ),
            stage + " method " + method
        )
    for i, comparator in enumerate(METHODS[1:]):
        x = effects(stage, "GAIT_CALIBRATED", comparator)
        require(
            close(median(x), TABLE[stage]["effects"][i]),
            stage + " paired " + comparator
        )
        require(
            int((x > 0).sum()) ==
            TABLE[stage]["positive"][i],
            stage + " positive patient count"
        )

print("TABLE 2 METHODS/PAIRED EFFECTS: PASS")

# B. Twenty independent Wilcoxon/Holm comparisons.
record = load("multiplicity.csv")
require(len(record) == 20, "multiplicity row count")
families = defaultdict(list)

for row in record:
    x = effects(
        row["stage"],
        row["better_method"],
        row["comparator_method"]
    )
    p = signed_p(x)
    require(
        pclose(p, row["raw_p"]),
        "E7 raw p " + row["contrast"]
    )
    families[row["family"]].append((row, p))

for family, entries in families.items():
    expected_size = int(entries[0][0]["family_size"])
    require(
        len(entries) == expected_size,
        family + " size"
    )
    adjusted = holm([p for _, p in entries])
    for (row, _), actual in zip(entries, adjusted):
        require(
            pclose(actual, row["holm_adjusted_p"]),
            family + " Holm " + row["contrast"]
        )

print("E7 WILCOXON/HOLM: 20/20 PASS")

# C. Six Table 2 patient-bootstrap confidence intervals.
ci_count = 0

for row in record:
    if row["better_method"] != "GAIT_CALIBRATED":
        continue
    if row["comparator_method"] not in METHODS[1:]:
        continue
    if row["family"] not in ("original_e7", "primary_new"):
        continue

    x = effects(
        row["stage"],
        row["better_method"],
        row["comparator_method"]
    )
    lo, hi = bootstrap(x, 20260829)
    require(
        close(lo, row["ci_low_mm"], 0.005)
        and close(hi, row["ci_high_mm"], 0.005),
        "E7 bootstrap " + row["contrast"]
    )
    ci_count += 1

require(ci_count == 6, "six Table 2 confidence intervals")
print("TABLE 2 BOOTSTRAP 95%CI: 6/6 PASS")

# D. E6 main perturbation point estimates.
perturbations = [
    ("position.csv", "amplitude_mm", 10.0,
     {"e1": 4.69, "e4": 0.61}),
    ("orientation.csv", "amplitude_deg", 5.0,
     {"e1": 19.83, "e4": 1.87}),
]

for filename, field, amplitude, expected in perturbations:
    rows = load(filename)
    for stage, target in expected.items():
        selected = [
            float(r["median_change_mm"])
            for r in rows
            if r["stage"].lower() == stage
            and "BIAS" in r["noise_model"].upper()
            and abs(float(r[field])-amplitude) < 1e-8
        ]
        require(len(selected) == 50, filename + " cohort")
        require(
            close(median(selected), target),
            filename + " " + stage
        )

# E6 orientation repetition heterogeneity.
repetitions = load("orientation_repeats.csv")
by_patient = defaultdict(list)
for r in repetitions:
    if (
        r["stage"].lower() == "e1"
        and "BIAS" in r["noise_model"].upper()
        and abs(float(r["amplitude_deg"])-5) < 1e-8
    ):
        by_patient[r["patient"]].append(
            float(r["median_anthro_advantage_mm"])
        )

require(len(by_patient) == 50, "E6 repeated patients")
require(
    all(len(x) == 10 for x in by_patient.values()),
    "E6 10 repeats"
)
always = sum(all(v > 0 for v in x)
             for x in by_patient.values())
never = sum(all(v <= 0 for v in x)
            for x in by_patient.values())
mixed = 50-always-never
require(
    (always, never, mixed) == (14, 13, 23),
    "E6 heterogeneity"
)
print("E6 PERTURBATION/REPETITION CHECK: PASS")

# E. E6 matched temporal comparisons.
def fraction_value(text):
    x = float(text)
    return x if x <= 1 else x/100

def mode_name(text):
    text = text.upper()
    if "RANDOM" in text:
        return "RANDOM"
    if "WINDOW" in text or "CONTIGUOUS" in text:
        return "WINDOW"
    raise RuntimeError("Unknown sampling mode")

paired = defaultdict(dict)

for r in load("information_repeats.csv"):
    f = fraction_value(r["fraction"])
    key = (
        round(f, 4),
        r["endpoint"].lower(),
        r["patient"],
        r["repeat"]
    )
    m = mode_name(r["sampling_mode"])
    require(m not in paired[key], "duplicate sampling")
    paired[key][m] = float(r["median_delta_mm"])

patient_effects = defaultdict(
    lambda: defaultdict(list)
)

for (fraction, endpoint, patient, repeat), values in paired.items():
    require(
        set(values) == {"WINDOW", "RANDOM"},
        "missing matched sampling mode"
    )
    patient_effects[
        (fraction, endpoint)
    ][patient].append(
        values["WINDOW"]-values["RANDOM"]
    )

frozen_e6 = load("information_paired.csv")
require(len(frozen_e6) == 8, "eight E6 comparisons")
e6_records = []

for row in frozen_e6:
    f = fraction_value(row["fraction"])
    endpoint = row["endpoint"].lower()
    grouped = patient_effects[(round(f, 4), endpoint)]
    require(len(grouped) == 50, "E6 patient count")
    require(
        all(len(x) == 10 for x in grouped.values()),
        "E6 ten repeats per patient"
    )

    # Pandas original groupby sorts patient IDs.
    values = np.asarray([
        median(grouped[p])
        for p in sorted(grouped)
    ])

    actual_effect = median(values)
    actual_p = signed_p(values)

    require(
        close(
            actual_effect,
            row["patient_median_difference_mm"],
            0.001
        ),
        "E6 matched patient median"
    )
    require(
        pclose(actual_p, row["wilcoxon_raw_p"]),
        "E6 Wilcoxon"
    )

    seed = (
        20261024 + int(f*1000)
        + (0 if endpoint == "e1" else 10000)
    )
    lo, hi = bootstrap(values, seed)

    require(
        close(lo, row["bootstrap95_low_mm"], 0.005)
        and close(hi, row["bootstrap95_high_mm"], 0.005),
        "E6 bootstrap"
    )
    e6_records.append((row, actual_p))

corrected = holm([p for _, p in e6_records])

for (row, _), value in zip(e6_records, corrected):
    require(
        pclose(value, row["holm8_p"]),
        "E6 Holm8"
    )

print("E6 TEMPORAL EFFECT/CI/HOLM: 8/8 PASS")
print("ALL STATISTICAL CHECKS PASSED")
