import argparse
import pandas as pd
import matplotlib.pyplot as plt
import ast
import numpy as np

def parse_data(input_file):
    try:
        # Read data from the input file
        data = pd.read_csv(input_file, sep=":", header=None, names=["X", "Molecules"])
    except FileNotFoundError:
        print("File not found. Please provide a valid input file.")
        return
    except pd.errors.ParserError:
        print("Invalid file format. Please ensure the input file is in the correct format.")
        return

    # Use ast.literal_eval to parse the nested list
    data["Molecules"] = data["Molecules"].apply(lambda x: ast.literal_eval(x))

    return data

def plot_data(input_file, output_file):
    data = parse_data(input_file)

    x_values = [1, 2, 3, 4, 5]  # Updated x-axis values
    energy = [5, 5.5, 6, 6.5, 7]  # Updated x-axis values
    y_values = data["Molecules"]

    # Create a list of colors for the box plots
    colors = ['#3D9970', '#FF4136', '#FF851B', '#FFC300', '#FF85FF']

    plt.figure(figsize=(12, 6))

    box_width = 0.15  # Adjust the width of the boxes

    for i, sublist in enumerate(y_values):
        temp = [-0.4, -0.2, 0, 0.2, 0.4]  # Adjust the position for each box plot
        positions = [item + i + 1 for item in temp]
        labels = [str(energy[i])] * 5  # Each x-axis value has 5 boxes

        boxprops = {'edgecolor': 'k'}  # Box outline color

        # Use the updated boxprops parameter
        boxplot = plt.boxplot(sublist, positions=positions, labels=labels, patch_artist=True, boxprops=boxprops)

        # Set the facecolor for each box in the current sublist
        for box, color in zip(boxplot['boxes'], colors):
            box.set(facecolor=color)

    plt.xlabel("Quantile")
    plt.ylabel("Explained Intensities")
    plt.title("Explained Intensities of 5 Ecom vs. Quantile of Predicted Intensities")
    plt.grid(True)

    # Create legend
    box_plots = [plt.Rectangle((0, 0), 1, 1, fc=colors[i]) for i in range(len(y_values))]    
    plt.legend(box_plots, energy)

    # Set the x-axis ticks to the updated values
    plt.xticks(np.arange(1, 6), x_values)

    # Save the plot as an image file
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot box plots from a file")
    parser.add_argument("input_file", type=str, help="Path to the input file")
    parser.add_argument("output_file", type=str, help="Path to the output image file")

    args = parser.parse_args()
    plot_data(args.input_file, args.output_file)
