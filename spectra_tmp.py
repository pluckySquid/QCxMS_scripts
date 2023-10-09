import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# Define the data for the heatmap
molecules = ['Molecule A', 'Molecule B', 'Molecule C', 'Molecule D']
# Create a random 7x18 list
data = np.random.rand(7, 18)
print(data)

# Create heatmap using seaborn
sns.heatmap(data, cmap='YlGnBu')

plt.title('Intensity of molecules')
plt.xlabel('Molecules')
plt.ylabel('Intensity level')

# Display the plot
plt.savefig('spec.png')
plt.show()
