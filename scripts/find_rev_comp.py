from Bio import SeqIO
from Bio.Seq import Seq

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def count_complements(fastq_file_1, fastq_file_2, allowed_errors=2):
    # FASTQ file sets #
    sequences_1 = []
    sequences_2 = set()  # Use a set for faster searching #

    # Read the first FASTQ file #
    for record in SeqIO.parse(fastq_file_1, "fastq"):
        sequences_1.append(str(record.seq))

    # Read the second FASTQ file #
    for record in SeqIO.parse(fastq_file_2, "fastq"):
        sequences_2.add(str(record.seq))

    # Initialize counters #
    normal_matches = 0
    reverse_complement_matches = 0

    # Check for normal and reverse complements #
    for seq in sequences_1:
        reverse_complement = str(Seq(seq).reverse_complement())
        
        for target_seq in sequences_2:
            # Only compare sequences of equal length #
            if len(seq) == len(target_seq):
                if hamming_distance(seq, target_seq) <= allowed_errors:
                    normal_matches += 1
                    break  # Avoid double-counting this match #
                if hamming_distance(reverse_complement, target_seq) <= allowed_errors:
                    reverse_complement_matches += 1
                    break

    return normal_matches, reverse_complement_matches

if __name__ == "__main__":
    fastq_file_2 = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/processed_reads/trimmed/A7480_S88_L001.merged.fq"
    fastq_file_1 = "/home/mikhail/Code/MFD-ILP/brendans_data_and_code/data/BaseCalls/A7480_S88_L001_R1_001.fastq"
    allowed_errors = 1  # You can change this to allow more errors
    normal, reverse_comp = count_complements(fastq_file_1, fastq_file_2, allowed_errors)
    print(f"Normal matches (with {allowed_errors} errors): {normal}")
    print(f"Reverse complement matches (with {allowed_errors} errors): {reverse_comp}")
