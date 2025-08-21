#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# === Read data from dissipation.txt ===
data = np.loadtxt("dissipation.txt")

# Extract columns
time = data[:, 0]
entropy = data[:, 1]
cost = data[:, 2]
potential = data[:, 3]

# === Save to CSV ===
df = pd.DataFrame({
    "time": time,
    "entropy": entropy,
    "cost": cost,
    "potential": potential
})
df.to_csv("dissipation.csv", index=False)

print("Saved dissipation.csv with columns: time, entropy, cost, potential")

# === Plot ===
plt.figure(figsize=(8, 5))
plt.plot(time, entropy, label="Entropy", lw=2)
plt.plot(time, cost, label="Kinetic Cost", lw=2)
plt.plot(time, potential, label="Dissipation Potential", lw=2)

plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Potentials over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the plot
plt.savefig("dissipation_plot.png", dpi=300)
print("Saved dissipation_plot.png")

# Show the plot
plt.show()

