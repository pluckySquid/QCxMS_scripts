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
    energy = [30, 40, 50, 60, 70, 80, 90, 100]  # Updated x-axis values
    y_values = data["Molecules"]
    print(y_values)

    # Create a list of unique colors for the box plots
    colors = ['#3D9970', '#FF4136', '#FF851B', '#FFC300', '#FF85FF', '#FF5733', '#3D86BC', '#FFAC45']

    plt.figure(figsize=(12, 6))

    bar_width = 0.1  # Adjust the width of the boxes

    for i, sublist in enumerate(y_values):
        temp = [-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35]  # Adjust the position for each box plot
        positions = [item + i + 1 for item in temp]
        labels = [30, 40, 50, 60, 70, 80, 90, 100] #[str(energy[i])]   # Each x-axis value has 8 boxes
        print("labels", labels)

        boxprops = {'edgecolor': 'k'}  # Box outline color

        print("sublist: ", len(sublist))
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
    box_plots = [plt.Rectangle((0, 0), 1, 1, fc=colors[i]) for i in range(len(energy))]
    plt.legend(box_plots, energy)

    # Set the x-axis ticks to the updated values
    plt.xticks(np.arange(1, 6), x_values)

    # Save the plot as an image file
    plt.savefig(output_file)
    plt.show()

def plot_data_normalized(input_file, output_file):
    output_file = "normalized_" + output_file 
    data = parse_data(input_file)

    x_values = [1, 2, 3, 4, 5]  # Updated x-axis values
    energy = [30, 40, 50, 60, 70, 80, 90, 100]  # Updated x-axis values
    y_values = data["Molecules"]
    first_quantile = y_values[0]

    y_values_normalized = []
    for quantile_list in y_values:
        temp_list_of_list = []
        for inner_list, first_energy_list in zip(quantile_list, first_quantile):
            divided_list = []
            for elem_a, elem_b in zip(inner_list, first_energy_list):
                if elem_b == 0:
                    divided_list.append(0)
                else:
                    divided_list.append(elem_a/elem_b)
            temp_list_of_list.append(divided_list)    
        y_values_normalized.append(temp_list_of_list)

    y_values = y_values_normalized
    print("y_values: ", y_values)
    # Create a list of unique colors for the box plots
    colors = ['#3D9970', '#FF4136', '#FF851B', '#FFC300', '#FF85FF', '#FF5733', '#3D86BC', '#FFAC45']

    plt.figure(figsize=(12, 6))

    bar_width = 0.1  # Adjust the width of the boxes

    for i, sublist in enumerate(y_values):
        temp = [-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35]  # Adjust the position for each box plot
        positions = [item + i + 1 for item in temp]
        labels = [30, 40, 50, 60, 70, 80, 90, 100] #[str(energy[i])]   # Each x-axis value has 8 boxes
        print("labels", labels)

        boxprops = {'edgecolor': 'k'}  # Box outline color

        print("sublist: ", len(sublist))
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
    box_plots = [plt.Rectangle((0, 0), 1, 1, fc=colors[i]) for i in range(len(energy))]
    plt.legend(box_plots, energy)

    # Set the x-axis ticks to the updated values
    plt.xticks(np.arange(1, 6), x_values)

    # Save the plot as an image file
    plt.savefig(output_file)
    plt.show()

def plot_data_normalized_one_energy(input_file, output_file):
    output_file = "one_energy_" + output_file 
    data = parse_data(input_file)

    x_values = [1, 2, 3, 4, 5]  # Updated x-axis values
    energy = [30, 40, 50, 60, 70, 80, 90, 100]  # Updated x-axis values
    
    ith_energy = 6

    energy = [energy[ith_energy]]
    y_values = data["Molecules"]
    ith_in_y_values = []
    for list_of_list in y_values:
        ith_in_y_values.append([list_of_list[ith_energy]])
    y_values = ith_in_y_values
    print("2 y_values: ", y_values)

    first_quantile = y_values[0]

    y_values_normalized = []
    for quantile_list in y_values:
        temp_list_of_list = []
        for inner_list, first_energy_list in zip(quantile_list, first_quantile):
            divided_list = []
            for elem_a, elem_b in zip(inner_list, first_energy_list):
                if elem_b == 0:
                    divided_list.append(0)
                else:
                    divided_list.append(elem_a/elem_b)
            temp_list_of_list.append(divided_list)    
        y_values_normalized.append(temp_list_of_list)

    y_values = y_values_normalized

    # Create a list of unique colors for the box plots
    colors = ['#3D9970', '#FF4136', '#FF851B', '#FFC300', '#FF85FF', '#FF5733', '#3D86BC', '#FFAC45']

    plt.figure(figsize=(12, 6))

    bar_width = 0.1  # Adjust the width of the boxes

    for i, sublist in enumerate(y_values):
        temp = [0]  # Adjust the position for each box plot
        positions = [item + i + 1 for item in temp]
        labels = [30, 40, 50, 60, 70, 80, 90, 100] #[str(energy[i])]   # Each x-axis value has 8 boxes
        labels = [labels[ith_energy]]
        print("labels", labels)

        boxprops = {'edgecolor': 'k'}  # Box outline color

        print("sublist: ", len(sublist))
        # Use the updated boxprops parameter
        boxplot = plt.boxplot(sublist, positions=positions, labels=labels, patch_artist=True, boxprops=boxprops)

        # Set the facecolor for each box in the current sublist
        for box, color in zip(boxplot['boxes'], colors):
            box.set(facecolor=color)

    plt.xlabel("Quantile")
    plt.ylabel("Normalized Explained Intensities")
    plt.title("Normalized Explained Intensities by the first molecule VS. Quantile of Predicted Intensities")
    plt.grid(True)

    # Create legend
    box_plots = [plt.Rectangle((0, 0), 1, 1, fc=colors[i]) for i in range(len(energy))]
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
    plot_data_normalized(args.input_file, args.output_file)
    plot_data_normalized_one_energy(args.input_file, args.output_file)
