"""
Takes FASTQ files as input

•Constructs a de Bruijn graph from k-mers (with user-defined k)

•Identifies contigs by finding Eulerian paths

•Outputs contigs as a FASTA file
"""

"""
resources used:
https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf
"""
import argparse
from collections import defaultdict
from Bio import SeqIO


def read_fastq(file_path: str):
    """
    Read a FASTQ file and extract the sequences.
    
    Args:
        file_path: Path to the FASTQ file
        
    Returns:
        A list of sequences from the FASTQ file
    """
    reads = []
    with open(file_path, 'r') as f:
        for record in SeqIO.parse(f, 'fastq'):
            reads.append(str(record.seq))
    return reads


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process k-mers from a file')
    parser.add_argument('-k', '--kmer-size', type=int, required=True,
                      help='Size of k-mer')
    parser.add_argument('-f', '--file', type=str, required=True,
                      help='Path to input file')
    parser.add_argument('-o', '--output', type=str, default='contigs.fasta.dbg',
                      help='Path to output FASTA file')
    
    # Parse arguments
    args = parser.parse_args()
    data = read_fastq(args.file)
    
    # Call the processing function
    left_kmers, right_kmers = process_kmers(args.kmer_size, data)

    graph = create_de_bruijn_graph(left_kmers, right_kmers)

    contigs = find_eulerian_path(graph)

    write_contigs_to_fasta(contigs, args.output)
    print(f"Wrote {len(contigs)} contigs to {args.output}")


def find_eulerian_path(graph):
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    for u in graph:
        for v in graph[u]:
            out_deg[u] += 1
            in_deg[v] += 1

    start = None
    for v in set(in_deg) | set(out_deg):
        if out_deg[v] - in_deg[v] == 1:
            start = v
            break
    if not start:
        start = next(iter(graph))  # arbitrary start

    path, stack = [], [start]
    while stack:
        v = stack[-1]
        if graph[v]:
            stack.append(graph[v].pop())
        else:
            path.append(stack.pop())

    tour = path[::-1]
    return [tour[0] + ''.join(map(lambda x: x[-1], tour[1:]))]

    
def create_de_bruijn_graph(left_kmers: list[str], right_kmers: list[str]):
    graph = defaultdict(list)
    for left, right in zip(left_kmers, right_kmers):
        graph[left].append(right)
    return graph


def get_number_of_edges(graph: dict, left_node: str, right_node: str):
    return graph[left_node].count(right_node) if right_node in graph[left_node] else 0


def get_kmer_of_size(k: int, data: str):
    return [data[i:i+k] for i in range(len(data) - k + 1)]


def process_kmers(k: int, data: list[str]):
    """
    Process k-mers from the given file
    Args:
        k: Size of k-mer
        filepath: Path to input file
    """
    # Flatten the list of kmers
    kmers = [kmer for read in data for kmer in get_kmer_of_size(k, read)]
    left_kmers = []
    right_kmers = []
    for kmer in kmers:
        left_kmers.append(kmer[:-1])
        right_kmers.append(kmer[1:])
    return left_kmers, right_kmers


def write_contigs_to_fasta(contigs, output_file):
    """
    Write contigs to a FASTA file.
    
    Args:
        contigs: List of contig sequences
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig{i+1} length={len(contig)}\n")
            # Write sequence with 60 characters per line
            for j in range(0, len(contig), 60):
                f.write(f"{contig[j:j+60]}\n")

if __name__ == "__main__":
    main()
