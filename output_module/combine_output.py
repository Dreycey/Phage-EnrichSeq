import sys
import os
import argparse
import csv
import matplotlib.pyplot as plt
current_path = os.path.dirname(os.path.abspath(__file__))

'''
TODO: Get output from:
1. genome_comparison/clusters.csv
2. genome_comparison/abundances.csv
3. merge_overlap_filter/merge_overlap_out_refined.csv
4. merge_overlap_filter/merge_overlap_out.csv

'''

def plot_piechart(labels: list, abundances: list, filename: str):
    '''
    Plots a pie chart of relative abundances

    INPUT:
        list of "labels" (taxids or clusters)
        list of abundances
        file name

    OUTPUT:
        None. Displays a pie chart.
    '''
    # display phage names in legend?
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:

    fig1, ax1 = plt.subplots()
    ax1.pie(abundances, labels=labels, autopct='%1.2f%%')
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.savefig(filename)

def plot_barchart():
    '''
    Plots a bar chart of relative abundances

    INPUT:
    
    OUTPUT:
        None. Displays a bar chart.
    '''


def parse_csv(file: str):
    '''
    CSV format is 2 columns in all cases

    INPUT:
        file to parse

    OUTPUT:
        list of first column values
        list of abundances (second column)
    '''

    #TODO: Exception handling for file exists

    labels = []
    values = []
    try:
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                labels.append(row[0])
                values.append(float(row[1]))
    except FileNotFoundError:
        print(f'File {file} not found.')

    return labels, values


def parse_args(argv=None) -> argparse.Namespace:
    '''
    DESCRIPTION:
        This method takes in the arguments from the command and performs
        parsing.

    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-d", "--inputdir", help="parent directory of modules containing all CSVs to parse", required=True)
    parser.add_argument("-o", "--outputdir", help="desired location of all output files", required=True)
    return parser.parse_args(argv)


def main():
    # Parse abundances from CSV
    arguments = parse_args(argv=sys.argv[1:])
    clusters, cluster_abund = parse_csv(arguments.inputdir + '/genome_comparison/abundances.csv')

    # Display pie chart of cluster abundances
    plot_piechart(clusters, cluster_abund, arguments.outputdir + '/cluster_piechart.png')

    # Display bar chart of cluster abundances --> is this necessary?

    # TODO: can there be a single consolidated CSV file?



if __name__ == "__main__":
    main()

    