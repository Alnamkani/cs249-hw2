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
    parser.add_argument('-g', '--gfa', type=str, default=None,
                      help='Path to output GFA file (optional)')
    
    # Parse arguments
    args = parser.parse_args()
    data = read_fastq(args.file)
    
    # Call the processing function
    left_kmers, right_kmers = process_kmers(args.kmer_size, data)

    graph = create_de_bruijn_graph(left_kmers, right_kmers)

    # Export the graph to GFA if specified
    if args.gfa:
        export_graph_to_gfa(graph, args.gfa)
        print(f"Wrote graph to {args.gfa}")

    contigs = find_eulerian_path(graph)

    write_contigs_to_fasta(contigs, args.output)
    print(f"Wrote {len(contigs)} contigs to {args.output}")


def find_eulerian_path(graph):
    # Convert graph to a format that tracks edge multiplicity
    edge_count = defaultdict(lambda: defaultdict(int))
    for u in graph:
        for v in graph[u]:
            edge_count[u][v] += 1
    
    # Find connected components
    components = find_connected_components(graph)
    
    contigs = []
    
    # Process each connected component
    for component in components:
        if not component:
            continue
        
        # Calculate in and out degrees for each node
        in_deg = defaultdict(int)
        out_deg = defaultdict(int)
        for u in component:
            for v, count in edge_count[u].items():
                out_deg[u] += count
                in_deg[v] += count
        
        # Find a valid starting node
        start = None
        
        # First, try to find a node with out_deg > in_deg (Eulerian path start)
        for v in component:
            if out_deg[v] - in_deg[v] == 1:
                start = v
                break
        
        # If no such node exists, use any node with outgoing edges
        if not start:
            for v in component:
                if v in edge_count and any(edge_count[v].values()):
                    start = v
                    break
        
        # Skip if no valid start node found
        if not start:
            continue
        
        # Use Hierholzer's algorithm for finding Eulerian path
        path = []
        def dfs(node):
            for next_node in list(edge_count[node].keys()):
                if edge_count[node][next_node] > 0:
                    edge_count[node][next_node] -= 1
                    dfs(next_node)
            path.append(node)
        
        dfs(start)
        path.reverse()
        
        # Construct the contig from the path
        if len(path) > 1:
            contig = path[0] + ''.join(map(lambda x: x[-1], path[1:]))
            contigs.append(contig)
    
    # If no contigs were found, return the original behavior with one path
    if not contigs and graph:
        original_path = original_find_eulerian_path(graph)
        return original_path
    
    return contigs


def original_find_eulerian_path(graph):
    """Original implementation kept as fallback"""
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


def find_connected_components(graph):
    """Find all connected components in the graph"""
    # Create an undirected version of the graph for component analysis
    undirected = defaultdict(set)
    for u in graph:
        for v in graph[u]:
            undirected[u].add(v)
            undirected[v].add(u)
    
    visited = set()
    components = []
    
    # DFS to find connected components
    def dfs(node, component):
        visited.add(node)
        component.add(node)
        for neighbor in undirected[node]:
            if neighbor not in visited:
                dfs(neighbor, component)
    
    # Find all components
    for node in undirected:
        if node not in visited:
            component = set()
            dfs(node, component)
            components.append(component)
    
    return components


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


def export_graph_to_gfa(graph: dict, output_file: str):
    """Export a de Bruijn graph (dict {node:[neighbors]}) to GFA 1.0."""
    with open(output_file, "w") as f:
        f.write("H\tVN:Z:1.0\n")
        nodes = set(graph)
        for adj in graph.values(): nodes.update(adj)
        for n in nodes: f.write(f"S\t{n}\t{n}\n")   # segments
        for src, dsts in graph.items():             # links (+ orientation, 0‑overlap)
            for dst in dsts: f.write(f"L\t{src}\t+\t{dst}\t+\t0M\n")



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
