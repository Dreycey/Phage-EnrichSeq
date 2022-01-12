import csv
import matplotlib.pyplot as plt


'''
TODO: Get output from:
1. genome_comparison/clusters.csv
2. genome_comparison/abundances.csv
3. merge_overlap_filter/merge_overlap_out_refined.csv
4. merge_overlap_filter/merge_overlap_out.csv

'''

def plot_piechart():
    '''
    Plots a pie chart of relative abundances

    INPUT:
        list of "classes" (taxids or clusters)
        list of abundances
    OUTPUT:
        None. Displays a pie chart.
    '''
    # display phage names in legend?

def plot_barchart():
    '''
    Plots a bar chart of relative abundances

    INPUT:
    
    OUTPUT:
        None. Displays a bar chart.
    '''


def main():
    # Parse abundances from CSV
    # Display pie chart of cluster abundances
    # Display bar chart of cluster abundances
    clusters = []
    cluster_abund = []
    plot_piechart(clusters, cluster_abund)

if __name__ == "__main__":
    main()