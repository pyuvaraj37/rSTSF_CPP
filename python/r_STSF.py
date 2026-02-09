from aeon.classification.interval_based import RSTSF
from aeon.datasets import load_italy_power_demand
import json
import numpy as np

X_train, y_train = load_italy_power_demand(split="TRAIN")
X_test, y_test = load_italy_power_demand(split="TEST")

clf = RSTSF(n_estimators=10, n_intervals=5, random_state=0)  

clf.fit(X_train, y_train)  


y_pred = clf.predict(X_test)  

print(y_pred)
cnt = 0
for i in range(len(y_pred)):
    if y_pred[i] == y_test[i]:
        cnt += 1

print("Accuracy:", cnt / len(y_pred))
accuracy = cnt / len(y_pred)

test_data = {
    "dataset_name": "ItalyPowerDemand",
    "X_test": X_test.tolist(),
    "y_test": y_test.tolist(),
    "y_pred": y_pred.tolist(),
    "test_accuracy": float(accuracy),
    "num_samples": len(X_test),
    "series_length": X_test.shape[1] if len(X_test.shape) > 1 else len(X_test[0]),
    "num_classes": len(np.unique(y_train))
}
    
test_filename = "test_data.json"
with open(test_filename, 'w') as f:
    json.dump(test_data, f, indent=2)

