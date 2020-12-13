# Pyxis2_classes_2020
The module argparse allows for the addition of command line arguments and parsing. Command-line interfaces which serves as interface for the input makes it possible for the users to easily change default settings in the code.
The argparse module is initialized by creating the ArgumentParser object. The following steps are adding more arguments to the object through the add_argument method.


parser = argparse.ArgumentParser(description='Pyxis detects clusters of regulated genes and returns the confidence')

parser.add_argument('-v', metavar='v', type=str, help="GFF file") #GFF file

parser.add_argument('-s', metavar='s', type=str, help="List of regulated genes") #file with upregulated genes

parser.add_argument('-i', metavar='i', type=str, nargs='?', default='ID=gene:', const='ID=gene:', help='Search string used to identify Genes in GFF file (default: ID=gene:)')

parser.add_argument('-k', metavar='k', type=str, nargs='?', default='Name=', const='Name=', help='Search string used to identify Gene name in GFF file (default: Name=)')

parser.add_argument('-e', metavar='p', type=int, nargs='?', default=0.01, const=0.01, help='Maximal P-value used to make statistical test more rigorous (default: 0.05)')

parser.add_argument('-p', metavar='p', type=int, nargs='?', default=0.001, const=0.001, help='Lower P-value used to make statistical test more rigorous (default: 0.05)')

parser.add_argument('-l', metavar='l', type=int, nargs='?', default=5, const=5, help='Minimum window to detect gene clusters (default: 5)')
args = parser.parse_args()

print("The search string for the GFF file is: ", args.i)

print("The database string for the GFF file is: ", args.k)

print("The cut-off p-value is: ", args.p)

print("The minimum window for the detection of gene clusters is: ", args.l)

GFF_file = args.v

files = args.s

idstring = args.i

sig_value = args.p

The first parameter in the method is created to implement flags into the command-line and are defined as boolean operators. This also makes the interpretation of the command-line easier. The metavar parameter in the method provides a different name for optional arguments in the help message and secondly the type argument indicates the type of input parameter the command-line requires. Lastly the help parameter in the method serves as navigation when any of the input files are faulty or when the command-line is not implemented correctly. The first two arguments both require input files and are mandatory. The first parameter requires the GFF3 file for the particular species and the second parameter requires the text file containing the list of gene names. The last three parameters are optional, and if no value or strings are given the default parameters will be implemented.

The third parameter is a default string specified for a type of GFF3 file implemented, and is also the first optional parameter in the command-line. The string “ID=gene” can be changed to, for example, “ID=SGD”, as may be required for GFF3 files representing annotations from different genomes.  The other parameter k contains the default string “Name=” which is always present in lines for genes which are protein coding.  If the user opts to remove string from the parameter dubious ORFs and pseudo-genes can also be retrieved.

The fifth parameter in the command-line is used to set maximal significance parameter as cut-off p-value for any significant cluster. The sixth parameter is the value which will indicate the most significant clusters from the least significant. Lastly, the seventh parameter is the step size at which the code detects clusters.

These parameters serve as input to the Pyxis2 program. The parse_args() method inspects the command-line and convert each argument to the appropriate type and the command-line arguments are then set to each of the various parameters in the program.

The command-line arguments which Pyxis2 uses after the import of the Argparse module is:

python Pyxis2.py -v #GFF3_file -s #list of regulated_genes -i #optional
