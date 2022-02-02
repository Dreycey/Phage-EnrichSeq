import sys
import os
import argparse

# in house packages
import consolidator
import plotter


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
    parser.add_argument("-i", "--inputdir", help="parent directory of modules containing all CSVs to parse", required=True)
    parser.add_argument("-o", "--outputdir", help="desired location of all output files", required=True)
    return parser.parse_args(argv)


def main():
    arguments = parse_args(argv=sys.argv[1:])

    # consolidator tasks
    list_of_files = [arguments.inputdir + '/merge_overlap_filter/merge_overlap_out.csv',
                     arguments.inputdir + '/genome_comparison/cluster_abundances.csv',
                     arguments.inputdir + '/genome_comparison/cluster_members.csv']
    
    consolidator.copy_files(list_of_files, arguments.outputdir)

    consolidator.rename_file(arguments.outputdir + '/merge_overlap_out.csv',
                             arguments.outputdir + '/taxid_abundances.csv')

    # Create charts for cluster abundances
    clusters, cluster_abund = plotter.parse_csv(arguments.outputdir + '/cluster_abundances.csv')
    plotter.plot_piechart(clusters, cluster_abund, arguments.outputdir + '/cluster_pie.png')
    plotter.plot_barchart(clusters, cluster_abund, arguments.outputdir + '/cluster_bar.png')

    # Create charts for phage abundances
    phages, phage_abund = plotter.parse_csv(arguments.outputdir + '/taxid_abundances.csv')
    plotter.plot_piechart(phages, phage_abund, arguments.outputdir + '/phage_pie.png')
    plotter.plot_barchart(phages, phage_abund, arguments.outputdir + '/phage_bar.png')



if __name__ == "__main__":
    main()