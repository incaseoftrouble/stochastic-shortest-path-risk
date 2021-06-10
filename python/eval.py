import json
import sys

if __name__ == "__main__":
    with open(sys.argv[1], mode="rt") as f:
        data = json.load(f)

    if sys.argv[2] == "figure":
        threshold_data = {
            threshold: {
                "vi": {
                },
                "lp": {}
            } for threshold in data["thresholds"]
        }

        for threshold, values in data["result"]["vi"].items():
            threshold_data[threshold]["var"] = values["var"]
        for threshold, times in data["time"]["vi"]["thresholds"].items():
            threshold_data[threshold]["vi"] = times["total"]
        for threshold, times in data["time"]["lp"]["thresholds"].items():
            threshold_data[threshold]["lp"] = times["total"]

        for threshold in sorted(threshold_data.keys()):
            values = threshold_data[threshold]
            print(f"{threshold},{values['var']:d},{values['vi']:.3f},{values['lp']:.3f}")
