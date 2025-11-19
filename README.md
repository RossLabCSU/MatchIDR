# MatchIDR
MatchIDR is a program designed to identify protein regions with an amino acid composition that most closely matches the composition of a query protein sequence.

MatchIDR calculates the percent composition for each of the 20 canonical amino acids from a user-defined query sequence. It then exhaustively scans a proteome or set of proteins to find the region with the best compositional match within each protein. While MatchIDR is expected to be especially useful for intrinsically disordered regions (IDRs) as query sequences, the mathematical process is identical and equally effective for non-IDR query sequences.

## Dependencies
MatchIDR requires NumPy version 1.20.0 or higher.

## Basic Usage
    python MatchIDR.py Proteome_File Query_Sequences_File [-o OUTPUT_FILE] [-n MINIMUM_WINDOW] [-x MAXIMUM_WINDOW] [-m METHOD] [-s]

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Proteome_File | The name of the file containing the protein sequences you wish to search for compositional matches (in FASTA format). |
| Query_Sequences_File | The name of the file containing the query sequence(s) in FASTA format (i.e., the reference sequences that you wish to find compositional matches for). |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -o OUTPUT_FILENAME | --output_file OUTPUT_FILENAME | Name of the output file that will contain the MatchIDR results (tab-separated values format). If no output filename is specified, the output file will be named based with a random job ID and the suffix "_MatchIDR_Results.tsv". |
| -n MINIMUM_WINDOW_SIZE | --minimum_window MINIMUM_WINDOW_SIZE | Smallest sliding window size to use during the MatchIDR proteome scan. |
| -x MAXIMUM_WINDOW_SIZE | --maximum_window MINIMUM_WINDOW_SIZE | Largest sliding window size to use during the MatchIDR proteome scan. |
| -m DISTANCE_METHOD | --method DISTANCE_METHOD | Distance method used to calculate compositional identity between two protein sequences. Manhattan distance (default) is recommended, as it makes the fewest assumptions of important protein features and enables consistent calculation of compositional identity from compositional distances. |
| -s | --seqset_average | Use the average composition vector for the set of sequences provided in the Query_Sequences_File for the distance calculations. This effectively represents your query sequences as a single "average" sequence. |
<br />

### Run MatchIDR with Default Parameters
MatchIDR is a Python3 script designed to be run as a stand-alone command-line application. To run MatchIDR on your sequences of interest, download the MatchIDR.py script and save to a location containing a FASTA file with your query sequence(s) and a FASTA file with your proteome of interest (or set of sequences). Navigate to that location via command line, and run MatchIDR with the following command (will use default parameters):

    python MatchIDR.py Proteome_File Query_Sequences_File

NOTE: Make sure to include the file extension in the command above for your file containing FASTA-formatted sequences. FASTA files will often have the file extension ".fa", ".fsa", or ".fasta", but are sometimes also provided as plain-text files (.txt), which should still work with MatchIDR. MatchIDR is designed to output your results in a **t**ab-**s**eparated **v**alues (.tsv) file.

## Detailed Usage and Customizable Parameters
The following sections illustrate the usage of optional command-line arguments.
<br/>
### Output filename (-o)
By default, MatchIDR will create an output file named with a random job ID and the suffix "_MatchIDR_Results.tsv". To specify an alternative output filename:

    python MatchIDR.py Proteome_File Query_Sequences_File -o Output_Filename

NOTE: MatchIDR will automatically add the ".tsv" extension to any output filename that you provide.

### Minimum window size (-n) and maximum window size (-x)
By default, MatchIDR uses a smallest sliding window size of 50 amino acids and a largest sliding window size of 200 amino acids. To use an alternative window size, use the ```-n``` and/or ```-x``` flags, with each flag followed by any positive integer value. For example, a search using a window size range of 80-120:

    python MatchIDR.py Proteome_File Query_Sequences_File -n 80 -x 120

### Distance method (-m)
The compositions of any two sequences are compared using a mathematical distance between the two corresponding composition vectors. By default, the distance method used by MatchIDR is the Manhattan distance (also known as the "city block" distance). Euclidean distance is also offered as an option for calculating distances:

    python MatchIDR.py Proteome_File Query_Sequences_File -m euclidean

### Sequence set average (-s)
Rather than using a single query sequence to find compositional matches in a proteome, there are cases where a user may want to perform a composition-matching search based on the average compositional characteristics of a set of sequences. For example, a user may have a set of IDR sequences that all localize to stress granules and have similar compositional features. Rather than performing individual searches with each IDR as a query sequence, the search can be performed by calculating the average compositional characteristics from the query sequence set ("Query_Sequences_File"), then using this as a representation of a single query sequence:

    python MatchIDR.py Proteome_File Query_Sequences_File --seqset_average

## License info
MatchIDR is subject to the terms of the MIT license.
