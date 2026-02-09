from aeon.classification.interval_based import RSTSF
from aeon.datasets import load_italy_power_demand
import json
import numpy as np
from aeon.utils.numba.stats import row_mean, row_median, row_std, row_slope, row_numba_min, row_numba_max, row_iqr, row_count_mean_crossing, row_count_above_mean

func_map = {
    row_mean: 0,
    row_median: 1,
    row_std: 2,
    row_slope: 3,
    row_numba_min: 4,
    row_numba_max: 5,
    row_iqr: 6,
    row_count_mean_crossing: 7,
    row_count_above_mean: 8,
}

X_train, y_train = load_italy_power_demand(split="TRAIN")
X_test, y_test = load_italy_power_demand(split="TEST")

clf = RSTSF(n_estimators=10, n_intervals=5, random_state=0)  

clf.fit(X_train, y_train)  


X_Diff = clf._series_transformers[0].transform(X_test)
X_Per = clf._series_transformers[1].transform(X_test)
X_Ar = clf._series_transformers[2].transform(X_test)

all_caf = []

for si in clf._transformers:
    all_caf.append([])
    for interval in si.intervals_:
        copy_interval = list(interval)
        copy_interval[3] = func_map.get(copy_interval[3], 99)  # Map function to its corresponding number
        all_caf[-1].append(copy_interval)

def py(x):
    return x.item() if hasattr(x, "item") else x

trees = []
for tree_idx, tree in enumerate(clf.clf_.estimators_):
    tree_ = tree.tree_
    trees.append([])
    for i in range(tree_.node_count):
        feature = tree_.feature[i]
        threshold = tree_.threshold[i]
        left = tree_.children_left[i]
        right = tree_.children_right[i]
        values = tree_.value[i][0]  # Shape is (1, n_classes)
        #print(f"Node {i}: feature={feature}, threshold={threshold}, left={left}, right={right}, values={values}")
        trees[-1].append([
            py(feature),
            py(threshold),
            py(left),
            py(right),
            [py(v) for v in values]
        ])

y_pred = clf.predict(X_test)  

cnt = 0
for i in range(len(y_pred)):
    if y_pred[i] == y_test[i]:
        cnt += 1

print("Accuracy:", cnt / len(y_pred))
accuracy = cnt / len(y_pred)

test_data = {
    "dataset_name": "ItalyPowerDemand",
    "X_test": X_test.squeeze(axis=1).tolist(),
    "X_Diff": X_Diff.squeeze(axis=1).tolist(),
    "X_Per": X_Per.squeeze(axis=1).tolist(),
    "X_Ar": X_Ar.squeeze(axis=1).tolist(),
    "y_test": y_test.astype(int).tolist(),
    "y_pred": y_pred.astype(int).tolist(),
    "all_candidate_agg_feats": all_caf.squeeze(axis=1).tolist(),
    "trees": trees,
    "test_accuracy": float(accuracy),
    "num_samples": len(X_test),
    "series_length": X_test.shape[1] if len(X_test.shape) > 1 else len(X_test[0]),
    "num_classes": len(np.unique(y_train))
}
    
test_filename = "test_data.json"
with open(test_filename, 'w') as f:
    json.dump(test_data, f, indent=2)

