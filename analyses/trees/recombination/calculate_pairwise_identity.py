import argparse
from Bio import SeqIO

def calculate_percent_identity(seq1, seq2, window_size):
    """Calculates the percent identity between two sequences in a sliding window."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length for pairwise comparison.")
    
    num_windows = len(seq1) - window_size + 1
    result = []

    # Slide the window across the sequences
    for start in range(num_windows):
        window_seq1 = seq1[start:start + window_size]
        window_seq2 = seq2[start:start + window_size]

        # Count matches
        matches = sum(1 for a, b in zip(window_seq1, window_seq2) if a == b)
        percent_identity = (matches / window_size) * 100

        # Store the result: window start position and percent identity
        result.append((start + 1, percent_identity))  # start + 1 for 1-based indexing

    return result

def process_fasta_alignment(file_path, window_size):
    """Process the FASTA file and calculate pairwise percent identity."""
    with open(file_path, 'r') as file:
        # Read sequences from the FASTA file
        sequences = list(SeqIO.parse(file, "fasta"))
        
        if len(sequences) != 2:
            raise ValueError("The input FASTA file must contain exactly two sequences.")
        
        seq1 = str(sequences[0].seq)
        seq2 = str(sequences[1].seq)
        
        # Calculate percent identity in sliding windows
        return calculate_percent_identity(seq1, seq2, window_size)

def write_output(results, output_file):
    """Writes the results to a tab-delimited file."""
    with open(output_file, 'w') as file:
        file.write("Window Start\tPercent Identity\n")
        for start, percent_identity in results:
            file.write(f"{start}\t{percent_identity:.2f}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate pairwise percent identity between two nucleotide sequences in sliding windows.")
    parser.add_argument("input_fasta", type=str, help="Input FASTA file containing two sequences for alignment.")
    parser.add_argument("window_size", type=int, help="The sliding window size for calculating percent identity.")
    parser.add_argument("output_file", type=str, help="Output file to store the percent identity results.")

    args = parser.parse_args()

    # Process the alignment and calculate percent identity
    results = process_fasta_alignment(args.input_fasta, args.window_size)
    
    # Write the results to the specified output file
    write_output(results, args.output_file)
    
    print(f"Output written to {args.output_file}")

if __name__ == "__main__":
    main()