import argparse
import requests
import urllib
import re
import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import textwrap

def plot_percentage(pathways, percentages):
    # Create a figure and axis object
    fig, ax = plt.subplots()

    # Set the title and axis labels
    ax.set_title('Percentage of Pathway Scores Above 0.6')
    ax.set_xlabel('Pathways')
    ax.set_ylabel('Percentage')

    # Create the bar chart
    ax.bar(pathways, percentages)

    # Add text labels for each bar
    for i, v in enumerate(percentages):
        ax.text(i, v+0.5, str(v) + '%', ha='center')

    # Save the plot to a file
    plt.savefig('pathway_scores_above_point_six.png')

    # Display the plot
    plt.show()

def plot_heatmap(pathway_scores_dict):
    # Determine the maximum length of the score lists
    max_length = max(len(scores) for scores in pathway_scores_dict.values())

    # Create a 2D numpy array filled with NaN values
    heatmap_data = np.empty((max_length, len(pathway_scores_dict)))
    heatmap_data[:] = np.nan

    # Fill the 2D numpy array with the values from the score lists
    for i, (molecule, scores) in enumerate(pathway_scores_dict.items()):
        scores = scores[::-1]  # Reverse the order of scores
        heatmap_data[:len(scores), i] = scores

    # Set up the figure and subplot configuration
    fig, ax = plt.subplots(figsize=(10, max_length))

    # Define a new colormap with inverted color mapping
    cmap = plt.cm.hot  # Choose the desired colormap
    new_cmap = colors.ListedColormap(cmap(np.flipud(np.linspace(0, 1, 256))))

    # Create the heatmap plot with the new colormap
    im = ax.imshow(heatmap_data, cmap=new_cmap, interpolation='nearest', aspect='auto')

    # Invert the y-axis
    ax.invert_yaxis()

    # Configure the plot axes
    ax.set_xticks(range(len(pathway_scores_dict)))
    ax.set_xticklabels([])  # Remove x-axis tick labels
    ax.set_xlabel('')

    ax.set_yticks(range(max_length))
    ax.set_yticklabels(range(1, max_length + 1))  # Set y-axis tick labels as score indices
    ax.set_ylabel('Score Index')
    # Change the font size of the y-axis labels
    ax.yaxis.get_label().set_size(15)

    # Add a colorbar
    cbar = fig.colorbar(im)
    cbar.set_label('Score')

    # Modify the font properties of the colorbar
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=25, weight='bold')
    cbar.ax.set_ylabel('Score', fontsize=25, weight='bold')

    # Adjust the spacing between subplots
    fig.tight_layout()

    # Add molecules outside the plot with line wrapping
    molecules = list(pathway_scores_dict.keys())
    x_positions = np.arange(len(molecules))
    y_position = -0.7
    for x, molecule in zip(x_positions, molecules):
        wrapped_lines = textwrap.wrap(molecule, width=10)  # Adjust the width as needed
        displayed_text = '\n'.join(wrapped_lines)
        ax.text(x, y_position, displayed_text, transform=ax.transData,
                ha='center', va='top', fontsize=15)

    # Save the plot as an image file (e.g., PNG or PDF)
    plt.savefig('heatmap.png', dpi=300, bbox_inches='tight')  # Change the filename and extension as desired

    # Display a message indicating successful saving
    print("Heatmap plot saved as 'heatmap.png'")
    plt.clf()

def main():

    molecules = []
    mol_score_dict = {}
    with open('results.csv') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            molecules.append(row[0].split('\t')[0])
            score = 0
            try:
                score = row[0].split('\t')[1]
            except:
                score = 0
            mol_score_dict[row[0].split('\t')[0]] = score
    print("mol_score_list: ", mol_score_dict)

    pathway_dict = {}
    for mol in mol_score_dict:
        smiles = urllib.parse.quote(mol)
        np_url = "https://npclassifier.ucsd.edu/classify?smiles={}".format(smiles) # NOTE: YOU need to URL encode the parameters or else sometimes it will fail
        r = requests.get(np_url)
    
        text = r.text

        match = re.search(r'"superclass_results": (\[.*?\])', text)

        if match:
            pathway_results = eval(match.group(1))
            print(len(pathway_results), mol, pathway_results)
            if len(pathway_results) == 0:
                if 'None' in pathway_dict:
                    pathway_dict['None'].append((mol, mol_score_dict[mol]))
                else:
                    pathway_dict['None'] = [(mol, mol_score_dict[mol])]
            else:
                for i in pathway_results:
                    if i in pathway_dict:
                        pathway_dict[i].append((mol, mol_score_dict[mol]))
                    else:
                        pathway_dict[i] = [(mol, mol_score_dict[mol])]
        else:
            print('No match found.')

    pathway_scores_dict = {}
    for pathway in pathway_dict:
        mol_list = pathway_dict[pathway]
        scores = [x[1] for x in mol_list]
        pathway_scores_dict[pathway] = scores
        #print(pathway, scores)
    print("pathway_scores_dict: ", pathway_scores_dict)

    plot_heatmap(pathway_scores_dict)

    pathway_scores_above_point_six_dict = {}
    for pathway in pathway_dict:
        scores = pathway_scores_dict[pathway]
        scores_above_point_six = [float(score) for score in scores if float(score) > 0.6]
        percentage_above_point_six = round((len(scores_above_point_six) / len(scores)) * 100, 2)
        pathway_scores_above_point_six_dict[pathway] = percentage_above_point_six
    print("pathway_scores_above_point_six_dict: ", pathway_scores_above_point_six_dict)
    
    plot_percentage(list(pathway_scores_above_point_six_dict.keys()), list(pathway_scores_above_point_six_dict.values()))

    # Create a list of labels and a list of values for each group
    labels = pathway_dict.keys()
    values = [[float(v) for v in pathway_scores_dict[k] if v] for k in pathway_scores_dict]

    # Create a figure and axis object
    fig, ax = plt.subplots()

    # Set the title and axis labels
    ax.set_title('Molecule Scores by np classification')
    ax.set_xlabel('class')
    ax.set_ylabel('Score')

    # Create the bar chart
    ax.bar(labels, [sum(v)/len(v) for v in values])

    # Add text labels for each bar
    for i, v in enumerate([sum(v)/len(v) for v in values]):
        ax.text(i, v+0.01, str(round(v,2)))

    # Save the plot to a file
    plt.savefig('classifier.png')

    # Display the plot
    plt.show()


if __name__ == "__main__":
    main()