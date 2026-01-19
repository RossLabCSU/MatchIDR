
"""
Description:
    MatchIDR is a composition-based method for identifying intrinsically disordered regions (IDRs) 
    of greatest compositional similarity to a query IDR. 
    
    Refer to (LINK TO PAPER WHEN PUBLISHED) for a complete description of MatchIDR.

============================================================================================
License_info:
    MatchIDR is subject to the terms of the MIT license. For a complete description of 
    license terms, please see the license at https://github.com/RossLabCSU/MatchIDR.
"""

__author__ = 'Sean M Cascarina'
__copyright__ = 'Copyright 2025'
__credits__ = ['Sean M Cascarina']
__license__ = 'MIT'
__version__ = '1.0'
__maintainer__ = 'Sean M Cascarina'
__email__ = 'Sean.Cascarina@colostate.edu'


def main(args):

    # RUN get_params() FUNCTION TO EXTRACT AND VALIDATE SEARCH PARAMETERS FROM COMMAND LINE ARGUMENTS
    query_ids, query_seqs, query_percs, min_window, max_window, amino_acids, method, fasta_file, output_file = get_params(args)

    random_id = random.randint(0, 10000000)
    print('Start time:', str(datetime.datetime.now()))

    # PREP OUTPUT FILE
    output = open(output_file, 'w')
    output_params_header(output, random_id, query_ids, query_seqs, min_window, max_window, method, fasta_file)
    output.write('\t'.join( ['Protein Description','UniProt ID (when applicable)','Query Protein ID','Best Fragment Match','Compositional Identity', 'Distance Score for Best Fragment (' + method + ' Distance)','Domain Boundaries for Best Fragment'] + amino_acids) + '\n')
    
    # CALCULATE PROTEOME SIZE AND PERCENTAGE MILESTONES FOR TRACKING CODE EXECUTION PROGRESS
    proteome_size = get_proteome_size(fasta_file)
    total_seqs = 0
    perc_range = [int(proteome_size*x/10) for x in range(1, 10)]

    # MUST BE A LIST OF NP ARRAYS SINCE THE SHAPE WOULD BE INCONSISTENT (i.e., AS NP ARRAY OF NP ARRAYS)
    amino_acids = [np.array([aa]) for aa in amino_acids]
    
    # CONTAINER TO STORE DATA FOR SORTING PRIOR TO OUTPUT
    stored_output = []

    # LOOP OVER EACH SEQUENCE IN THE FASTA FILE
    h = open(fasta_file)
    for prot_id, seq in fasta_parser(h):
    
        # SPECIAL PARSING FOR UNIPROT PROTEOME FILES
        if prot_id.count('|') == 2:
            junk, uniprot, junk = prot_id.split('|')
        else:
            uniprot = 'N/A'
            
        # SKIP PROTEINS THAT ARE SHORTER THAN THE MINIMUM WINDOW SIZE
        if len(seq) < min_window:
            total_seqs += 1
            continue
            
        #PRINT PROGRESS AT 10% MILESTONES
        if total_seqs in perc_range:
            print(str(round(total_seqs/proteome_size*100)) + '% Complete' + '\t' + str(datetime.datetime.now()))
        
        #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
        if seq.count('*') == 1 and seq[-1] == '*':
            seq = seq[:-1]
        elif seq.count('*') > 1:
            print('\nRuntime warning: protID \"' + id[:15] + '...\" contains multiple stop codons, which will slightly affect the amino acid composition and separation calculations. Stop codons are automatically removed by the program from the C-terminus of each sequence but internal stop codons are not removed. Consider removing extra stop codons before evaluating, removing these sequences from analyses entirely, or evaluating sequences as-is.\n')
            if seq[-1] == '*':
                seq = seq[:-1]
                
        # FOR EACH PROTEIN IN THE PROTEOME, LOOP OVER THE QUERY PROTEINS AND FIND THE BEST MATCHING FRAGMENT. THE MATCHING IS PERFORMED BY THE run_exhaustive_compsim_search() FUNCTION.
        for query_index, query_perc_arr in enumerate(query_percs):
            query_id = query_ids[query_index]
            best_score, best_frags, best_sim_arrays = run_exhaustive_compsim_search(query_perc_arr, seq, min_window, max_window, amino_acids, method)
            comp_identity = 100 - best_score/2
            boundaries = []
            for best_frag in best_frags:
                start = seq.index(best_frag) + 1
                end = start + len(best_frag)
                bounds = '('+str(start) + '-' + str(end) + ')'
                boundaries.append(bounds)
            
            # STORE DATA FOR SORTING AND WRITING TO OUTPUT
            output_data = [str(x) for x in [prot_id, uniprot, query_id, ';'.join(best_frags), comp_identity, best_score, ';'.join(boundaries)] + best_sim_arrays[0]]
            stored_output.append( (comp_identity, output_data) )
            
        total_seqs += 1
        
    h.close()
    
    # WRITE STORED DATA TO OUTPUT FILE
    stored_output.sort(reverse=True)
    for output_data in stored_output:
        output_line = output_data[1]
        output.write('\t'.join(output_line) + '\n')
    output.close()
    
    print('End time:', str(datetime.datetime.now()))
    
    
def run_exhaustive_compsim_search(query_perc_arr, prot_seq, min_window, max_window, amino_acids, method):
    """Exhaustively scans a protein using window sizes from the minimum window to the maximum window (or the length of the protein).
        For each window, this function passes the window and query seq to the calc_similarity_Distance() function to 
        calculate the distance score for that window. After evaluating all windows, the best-scoring fragment, its associated
        similarity score, and contributions of the individual amino acids (or amino acid groups) to the final similarity score
        are returned. All of these values correspond to a single, exhaustively-searched protein.
        
    Returns:
        best_score (float): compositional distance score of best-matching fragment(s)
        best_frags (list of strings): sequence(s) of best-matching fragment(s). All sequences tied for the single best score are included in the list
        best_sim_arrays (list of lists): list of compositional distance array(s) for best matching fragment(s). Each sublist contains the mathematical contribution of each amino acid or amino acid group to the final composition distance score
    """

    # INITIALIZE VARIABLES WITH DUMMY VALUES
    best_score = 10000  
    best_frags = []
    best_sim_arrays = []
    for win_size in range(min_window, min(max_window+1, len(prot_seq)+1)):
        
        # CONVERT SEQUENCE TO NUMPY ARRAY
        seq_rep = np.array(list(prot_seq))
        
        # BREAK PROTEIN SEQUENCE INTO FRAGMENTS USING CURRENT WINDOW SIZE
        frags = sliding_window_view(seq_rep, win_size)
        
        # COUNT EACH AMINO ACID FOR EACH FRAG
        # CALCULATE COUNT OF EACH AMINO ACID FOR EVERY FRAGMENT
        # WILL BE ARRAY OF SUBARRAYS (20, n_frags), WHERE EACH SUBARRAY CONTAINS THE COUNT FOR ONE OF THE AMINO ACIDS FOR ALL FRAGS
        # NOTE THAT THIS NEEDS TO BE RESTRUCTURED BEFORE IT REPRESENTS THE COUNT ARRAY FOR EACH FRAGMENT
        # .sum(1) CALCULATES THE SUM OF THE COUNTS FOR EACH ROW (EACH ROW REPRESENTS A DIFFERENT SEQUENCE FRAGMENT)
        all_counts = np.array([np.isin(frags, aa).sum(1) for aa in amino_acids], dtype=np.float32)

        # RESTRUCTURES THE ARRAY OF SUBARRAYS SUCH THAT EACH ROW NOW REPRESENTS THE COUNT COMPOSITION OF EACH UNIQUE SEQUENCE FRAG
        # THIS EFFECTIVELY CONVERTS THE COLUMNS TO ROWS
        counts_arrs = np.column_stack(all_counts)

        # CONVERT COUNTS TO PERCENTS
        all_percs = counts_arrs / win_size * 100

        # CALCULATE DISTANCES
        if method.upper() == 'MANHATTAN':
            sim_arrays = np.abs(all_percs - query_perc_arr)
        elif method.upper() == 'EUCLIDEAN':
            sim_arrays = np.square(all_percs - query_perc_arr)

        scores = sim_arrays.sum(1)
        if method.upper() == 'EUCLIDEAN':
            scores = np.sqrt(scores)
        
        # GET MINIMUM SCORE, INDEX OF MINIMUM SCORE, AND THE SIMILARITY ARRAY FOR MINIMUM SCORE
        smallest_loc = np.argmin(scores)
        smallest_score = np.min(scores)
        smallest_sim_array = sim_arrays[smallest_loc]

        # STORE IF THAT IS THE BEST OVERALL SCORE SO FAR
        if smallest_score < best_score:
            best_score = smallest_score
            smallest_locs = np.where(scores==smallest_score)[0]     # NEED TO INDEX THE RESULT AT [0] BECAUSE THE RESULT IS A TUPLE WITH THE DESIRED NP ARRAY OF POSITIONS AS THE FIRST ITEM IN THE TUPLE (e.g., (array([408, 410, 411, 412, 413], dtype=int64),) )
            best_frags = [prot_seq[smallest_loc:smallest_loc+win_size] for smallest_loc in smallest_locs]
            best_sim_arrays = [list(sim_arrays[smallest_loc]) for smallest_loc in smallest_locs]

        elif smallest_score == best_score or np.isclose(smallest_score, best_score):
            smallest_locs = np.where(scores==smallest_score)[0]     # NEED TO INDEX THE RESULT AT [0] BECAUSE THE RESULT IS A TUPLE WITH THE DESIRED NP ARRAY OF POSITIONS AS THE FIRST ITEM IN THE TUPLE (e.g., (array([408, 410, 411, 412, 413], dtype=int64),) )
            best_frags += [prot_seq[smallest_loc:smallest_loc+win_size] for smallest_loc in smallest_locs]
            best_sim_arrays += [list(sim_arrays[smallest_loc]) for smallest_loc in smallest_locs]

    return best_score, best_frags, best_sim_arrays


def get_proteome_size(fasta_file):
    """Tally the # of proteins in the user-defined FASTA file.
    
    Returns:
        Proteome size (int)
    """
    h = open(fasta_file)
    proteome_size = 0
    for line in h:
        if line.startswith('>'):
            proteome_size += 1
            
    return proteome_size
    
    
def fasta_parser(file):
    """Parse each instance in a FASTA formatted file into a gene id and a sequence.
    
    Yields:
        id, seq (both strings)
    """
    #INITIALIZES GENE AND SEQ FOR FIRST ITERATION
    gene, seq = '', []
    
    #LOOPING THROUGH FILE USING GENERATOR
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
            
        #YIELDS GENE AND SEQ IF THE NEXT ID IS REACHED
        if line.startswith('>'):
            if gene != '':
                yield (gene[1:], ''.join(seq))
            gene, seq = line, []
        else:
            seq.append(line)
            
    #YIELDS FINAL INSTANCE IN FASTA FILE
    if gene != '':
        yield (gene[1:], ''.join(seq))
    
    
def calc_composition_Distance(string1, string2, amino_acids, method):
    """Function that calculates the total composition distance score based on one of two metrics.
    Manhattan Distance - This score is calculated by:
    1) calculating the % composition for an amino acid in the query sequence and subject sequence.
    2) calculating the absolute value of the difference in % compositions calculated in step 1.
    3) repeating steps 1-2 for all amino acids or amino acid groups.
    4) summing the values from step 3.
    
    Euclidean Distance - This score is determined by:
    1) calculating the % composition for an amino acid in the query sequence and subject sequence.
    2) calculating the square of the difference in % compositions calculated in step 1.
    3) repeating steps 1-2 for all amino acids or amino acid groups.
    4) summing the values from step 3.
    5) taking the square-root of the value from step 4.
    
    Returns:
        distance_score (float): compositional distance score
        distance_arr (list): compositional distance array containing the mathematical contribution of each amino acid or amino acid group to the final composition distance score
    """
    if method.upper() == 'MANHATTAN':
        distance_arr = [ abs( string1.count(aa)/len(string1) - string2.count(aa)/len(string2) ) *100 for aa in amino_acids ]
        distance_score = sum( distance_arr )
    elif method.upper() == 'EUCLIDEAN':
        distance_arr = [ ( ( string1.count(aa)/len(string1) - string2.count(aa)/len(string2) ) *100)**2 for aa in amino_acids ]
        distance_score = math.sqrt( sum( distance_arr ) )

    return distance_score, distance_arr
    
    
def calc_seqset_average(query_seqs, amino_acids):
    """Calculate the average percent composition for each amino acid from a set of query sequences.
    
    Returns:
        mean_comps (list): list of average percent composition values for the set of query sequences, in the same order as amino_acids.
    """
    
    df = {}
    for aa in amino_acids:
        comps = []
        for seq in query_seqs:
            comp = seq.count(aa) / len(seq) *100
            comps.append(comp)
        df[aa] = comps
        
    mean_comps = np.array( [sum(df[aa]) / len(df[aa]) for aa in amino_acids] )
    
    return mean_comps
    
    
def get_params(args):
    """Gather, validate, and define user parameters."""
    
    #GET USER-SPECIFIED PARAMETERS=================
    fasta_file = args.fasta_file

    # VALIDATE QUERY SEQUENCE
    query_seq_file = args.query_file
    try:
        h = open(query_seq_file)
        query_ids = []
        query_seqs = []
        for prot_id, seq in fasta_parser(h):
            
            is_invalid = sum([1 for x in set(seq) if x in 'ACDEFGHIKLMNPQRSTVWY']) != len(set(seq))
            if is_invalid:
                print('\nError:\nYour query IDR sequence must be in FASTA format and the sequence can only contain the 20 canonical amino acids as well as O (pyrrolysine) or U (selenocysteine). No other symbols are permitted.\n')
                exit()
                
            query_ids.append(prot_id)
            query_seqs.append(seq)
    except:
        print('\nError:\nYour query IDR sequence file could not be found. Common causes for this error include (but are not limited to):\n1) The query seqeuence file is not in the same folder/directory as this script,\n2) the query sequence file is misspelled,\n3) there are spaces in the name of the query sequence file, which must be handled correctly in the command (typically surrounded in quotes, but this may be OS-specific).\n')
        exit()
        
    # VALIDATE MINIMUM WINDOW SIZE
    min_window = args.minimum_window
    try:
        min_window = int(min_window)
        if min_window < 10:
            print('\nError:\nYour minimum window size must be an integer >9.\n')
            exit()
    except:
        print('\nError:\nYour minimum window size must be an integer >9.\n')
        exit()
        
    # VALIDATE MAXIMUM WINDOW SIZE
    max_window = args.maximum_window
    try:
        max_window = int(max_window)
        if max_window < 10 or min_window > max_window:
            print('\nError:\nYour maximum window size must be an integer >9 and must be larger than your minimum window size (default minimum is 20aa).\n')
            exit()
    except:
        print('\nError:\nYour maximum window size must be an integer >9 and must be larger than your minimum window size (default minimum is 20aa).\n')
        exit()
        
        
    all_amino_acids = 'ACDEFGHIKLMNPQRSTVWYOU'
    amino_acids = list(all_amino_acids)

    query_percs = []
    for seq in query_seqs:
        percs = np.array([seq.count(aa) / len(seq) * 100 for aa in amino_acids])
        query_percs.append(percs)
        
    use_seqset_average = args.seqset_average
    if use_seqset_average:
        query_percs = [calc_seqset_average(query_seqs, amino_acids)]
        query_ids = ['Sequence Average from '+query_seq_file]

    # VALIDATE DISTANCE METRIC
    method = args.method
    try:
        method = method.upper().replace('"', '')
        options = ['MANHATTAN', 'EUCLIDEAN']
        test = options.index(method)    # WILL FAIL IF THE USER-DEFINED METHOD IS NOT "MANHATTAN" OR "EUCLIDEAN".
    except:
        print('\nError:\nThe only available options for the method parameter are "Manhattan" or "Euclidean".\n')
        exit()
        
    output_file = args.output_file
    # USE RANDOM JOB ID IF NO OUTPUT FILE IS SPECIFIED. OTHERWISE, ADD .tsv EXTENSION TO PROVIDED FILENAME.
    if not output_file:
        random_id = random.randint(0, 10000000)
        output_file = f'{random_id}_MatchIDR_Results.tsv'
    else:
        if not output_file.endswith('.tsv'):
            output_file += '.tsv'
            
    # CHECK IF OUTPUT FILE ALREADY EXISTS
    try:
        h = open(output_file)
        h.close()
        user_input = input('\nWARNING: An output file with your specified output filename already exists. Do you wish to overwrite this file? (y/n)  ')
        if user_input.upper() == 'N' or user_input.upper() == 'NO':
            print('\nExiting without running MatchIDR...\n')
            exit()
    except SystemExit:
        exit()
    except:
        pass

    return query_ids, query_seqs, query_percs, min_window, max_window, amino_acids, method, fasta_file, output_file


def output_params_header(output, random_id, query_ids, query_seqs, min_window, max_window, method, fasta_file):
    """Write runtime parameters to the output file before running analyses.
    
    Returns:
        output (file handle)
    """
    
    output.write('>*RUNTIME PARAMETERS*\n')
    output.write('>Job ID: ' + str(random_id) + '\n')
    output.write('>Query IDR Sequence(s): ' + str(query_seqs) + '\n')
    output.write('>FASTA ID(s) of Query IDR(s): ' + str(query_ids) + '\n')
    output.write('>FASTA File(s): ' + fasta_file + '\n')
    output.write('>Minimum Window Size: ' + str(min_window) + '\n')
    output.write('>Maximum Window Size: ' + str(max_window) + '\n')
    output.write('>Compositional Similarity Method: ' + method + '\n\n')

    return output


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Automated IDR similarity search.\n')
    
    parser.add_argument('fasta_file', help="""Your sequence file (in FASTA format).""")
    
    parser.add_argument('query_file', help="""The name of the FASTA file containing your query IDR sequence(s).""")
                        
    parser.add_argument('-o', '--output_file', type=str, default=None,
                        help="""The name of the resulting output file. If no output filename is provided, the output file is named with a randomly generated job ID.
                        """)
                        
    parser.add_argument('-n', '--minimum_window', type=int, default=50,
                        help="""Minimum window size used in the search. Must be a whole number >9.
                        
                        default = 50
                        """)
                        
    parser.add_argument('-x', '--maximum_window', type=int, default=200,
                        help="""Maximum window size used in the search. Must be a whole number >9 and must be greater than or equal to the minimum window size.
                        NOTE: The maximum window size can be greater than or equal to the largest protein in your proteome file, but the effective maximum window size used in the search will never exceed the size of the protein being analyzed. A large difference between minimum window size and maximum window size may result in extremely long computation times.
                        
                        default = 200
                        """)
                        
    parser.add_argument('-m', '--method', type=str, default='manhattan',
                        help="""Method used for quantification of compositional similarity.
                        Available options are "Manhattan" and "Euclidean", which will calculate the Manhattan distance or Euclidean distance, respectively, for the composition vectors associated with your query IDR sequence and a protein sequence from the proteome.
                        """)

    parser.add_argument('-s', '--seqset_average', action='store_true',
                        help="""Searches for sequences that best match the composition vector
                        representing the average composition of the query sequences. This is essentially
                        like creating a single "average" query sequence from a set of multiple sequences,
                        then searching for fragments that best match this average sequence.""")

    args = parser.parse_args(arguments)
    
    return args


if __name__ == '__main__':
    import sys, argparse, random, datetime, math
    import numpy as np
    from numpy.lib.stride_tricks import sliding_window_view
    args = get_args(sys.argv[1:])
    main(args)