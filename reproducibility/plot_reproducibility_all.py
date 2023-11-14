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
    print("data", data)
    
    x_values = data["X"]
    y_values = data["Molecules"]

    # Create a list of colors for the box plots
    colors = ['#3D9970', '#FF4136', '#FF851B', '#FFC300', '#FF85FF']

    plt.figure(figsize=(12, 6))

    for i, sublist in enumerate(y_values):
        positions = np.arange(len(sublist)) + 1 + i*0.2  # Adjust the position for each box plot
        labels = [x_values[i]] * len(sublist)

        plt.boxplot(sublist, positions=positions, labels=labels, patch_artist=True,
                    boxprops={'facecolor': colors[i]})

    plt.xlabel("Number before ':'")
    plt.ylabel("Values")
    plt.title("Box Plots of Molecules")
    plt.grid(True)

    # Create legend
    box_plots = [plt.Rectangle((0, 0), 1, 1, fc=colors[i]) for i in range(len(y_values))]
    plt.legend(box_plots, x_values)

    # Save the plot as an image file
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot box plots from a file.")
    parser.add_argument("input_file", type=str, help="Path to the input file")
    parser.add_argument("output_file", type=str, help="Path to the output image file")

    args = parser.parse_args()
    plot_data(args.input_file, args.output_file)
