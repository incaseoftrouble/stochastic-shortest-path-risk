import json
import sys
from collections import defaultdict

if __name__ == "__main__":
    with open(sys.argv[1], mode="rt") as f:
        data = json.load(f)

    if sys.argv[2] == "figure":
        vi_thresholds = set(data["vi"]["results"].keys())
        lp_thresholds = set(data["lp"]["results"].keys())
        if vi_thresholds != lp_thresholds:
            raise ValueError

        threshold_data = {threshold: {} for threshold in vi_thresholds}
        for threshold, vi_values in data["vi"]["results"].items():
            threshold_data[threshold]["vi"] = vi_values["time"]
            threshold_data[threshold]["var"] = vi_values["var"]
        for threshold, lp_values in data["lp"]["results"].items():
            threshold_data[threshold]["lp"] = lp_values["time"]
            if lp_values["var"] != threshold_data[threshold]["var"]:
                raise ValueError(f'{lp_values["var"]} vs {threshold_data[threshold]["var"]} for {threshold}')

        for threshold, values in sorted(threshold_data.items(), reverse=True):
            print(f"{float(threshold):.3f},{values['var']:d},{values['vi']:.3f},{values['lp']:.3f}")
