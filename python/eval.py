import json
import sys
from collections import defaultdict

if __name__ == "__main__":
    with open(sys.argv[1], mode="rt") as f:
        data = json.load(f)

    if "vi" in data:
        vi_thresholds = set(data["vi"]["results"].keys())
        thresholds = vi_thresholds
    if "lp" in data:
        lp_thresholds = set(data["lp"]["results"].keys())
        thresholds = lp_thresholds

    if "vi" in data and "lp" in data:
        if vi_thresholds != lp_thresholds:
            raise ValueError
    threshold_data = {threshold: {} for threshold in vi_thresholds}
    if "vi" in data:
        for threshold, vi_values in data["vi"]["results"].items():
            threshold_data[threshold]["vi"] = vi_values["time"]
            threshold_data[threshold]["var"] = vi_values["var"]
            threshold_data[threshold]["cvar"] = vi_values["cvar"]
    if "lp" in data:
        for threshold, lp_values in data["lp"]["results"].items():
            threshold_data[threshold]["lp"] = lp_values["time"]
            threshold_data[threshold]["cvar"] = lp_values["cvar"]
            if lp_values["var"] != threshold_data[threshold]["var"]:
                raise ValueError(f'{lp_values["var"]} vs {threshold_data[threshold]["var"]} for {threshold}')

    for threshold, values in sorted(threshold_data.items(), key=lambda x: float(x[0]), reverse=True):
        print(f"{float(threshold):.10f},{values['var']:d},{values['cvar']:.5f},{values.get('vi', 0.0):.3f},{values.get('lp', 0.0):.3f}")
