import matplotlib.pyplot as plt
import numpy as np

# Your data
best_energy_values = [40, 40, 50, 70, 50, 30, 70, 50, 40, 60, 60, 60, 40, 40, 50, 100, 30, 30, 90, 30, 50, 70, 30, 70, 90, 30, 50, 30, 30, 90, 40, 50, 30, 90, 60, 60, 80, 90, 40, 90, 80, 70, 90, 90, 40, 90, 60, 70, 90, 90, 80, 30, 90, 80]

# Create an x-axis range from 1 to the number of elements in best_energy_values
x_axis_range = np.arange(1)

# Create a boxplot
plt.clf()
plt.boxplot(best_energy_values)
plt.title('Intensity Trend (Boxplot)')
plt.xlabel('Molecule Index')
plt.ylabel('Energy Value')
plt.xticks(x_axis_range, x_axis_range)  # Set x-axis tick labels
plt.savefig("each_mole_boxplot.png")
plt.show()  # Show the plot
