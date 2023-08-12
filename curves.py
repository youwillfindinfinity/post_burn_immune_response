import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Sigmoid function
def sigmoid(x, a, b):
    return 1 / (1 + np.exp(-a * (x - b)))

# Michaelis-Menten curve
def michaelis_menten(x, km, vmax):
    return (vmax * x) / (km + x)

# Parameter values for subplots
sigmoid_params = [(1, 4), (2, 4), (1, 5), (2, 5)]
mm_params = [(-96, -20), (-97, -25), (-100, -30), (-96, -33)]

# Create the figure and subplots using Seaborn
sns.set(style="whitegrid")
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Plot sigmoid subplots
for a, b in sigmoid_params:
    x = np.linspace(-10, 10, 100)
    y = sigmoid(x, a, b)
    y = y - min(y)  # Shift the curve to start from 0
    label = f'a = {a}, b = {b}'
    sns.lineplot(x=x, y=y, label=label, ax=axs[0])

# Plot Michaelis-Menten subplots
for km, vmax in mm_params:
    x = np.arange(0, 1000000, 1)  # Adjust the x range to 0-100
    y = michaelis_menten(x, km, vmax)
    label = f'km = {km}, vmax = {vmax}'
    sns.lineplot(x=x, y=y, label=label, ax=axs[1])

    # Set custom x-axis tick labels
    axs[1].set_xticks(np.arange(0, 110, 10))
    axs[1].set_xticklabels(np.arange(0, 110, 10))

# Set titles and legends
axs[0].set_title('Sigmoid Function')
axs[0].set_ylabel('Probability of activation')
axs[0].set_xlabel('Concentration of cytokines')
axs[0].set_xlim(0, 10)
axs[0].legend()

axs[1].set_title('Michaelis-Menten Curve')
axs[1].set_ylabel('m2deathspan')
axs[1].set_xlabel('Time(h)')
axs[1].set_xlim(95, 100)
axs[1].legend()

# Add a big title at the top
fig.suptitle('Comparison of Sigmoid and Michaelis-Menten Curves')

# Adjust layout
plt.tight_layout()

# Save the plot as a PNG
plt.savefig('curves_plot.png', dpi=300)

# Show the plot
plt.show()
