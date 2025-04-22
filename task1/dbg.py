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
    # Create a deep copy of the graph to avoid modifying the original
    graph_copy = defaultdict(list)
    for u in graph:
        graph_copy[u] = graph[u].copy()
    
    # Calculate in and out degrees for each node
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    for u in graph_copy:
        for v in graph_copy[u]:
            out_deg[u] += 1
            in_deg[v] += 1
    
    all_nodes = set(in_deg) | set(out_deg)
    
    # Find all connected components in the graph
    components = find_connected_components(graph_copy)
    
    contigs = []
    
    # Process each connected component
    for component in components:
        if not component:
            continue
            
        # Find a valid starting node for each component
        start = None
        
        # First, try to find a node with out_deg > in_deg (Eulerian path start)
        for v in component:
            if out_deg[v] - in_deg[v] == 1:
                start = v
                break
                
        # If no such node exists, try to find a node with equal in/out degree (Eulerian cycle)
        if not start:
            for v in component:
                if v in graph_copy and graph_copy[v]:
                    start = v
                    break
        
        # Skip if no valid start node found
        if not start:
            continue
            
        # Find path starting from this node
        path, stack = [], [start]
        while stack:
            v = stack[-1]
            if v in graph_copy and graph_copy[v]:
                stack.append(graph_copy[v].pop())
            else:
                path.append(stack.pop())
        
        # Construct the contig from the path
        if path:
            tour = path[::-1]
            if len(tour) > 1:  # Only add if we have a real path
                contig = tour[0] + ''.join(map(lambda x: x[-1], tour[1:]))
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


def export_graph_to_gfa(graph, output_file):
    """
    Export the de Bruijn graph to GFA format.
    
    Args:
        graph: De Bruijn graph as a dictionary
        output_file: Path to output GFA file
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments (nodes)
        segments = set()
        for src in graph:
            segments.add(src)
            for dst in graph[src]:
                segments.add(dst)
                
        for i, segment in enumerate(segments):
            f.write(f"S\t{segment}\t{segment}\n")
        
        # Write links (edges)
        for src in graph:
            for dst in graph[src]:
                # In GFA, links use the orientation (+/-) for each segment
                # For simplicity, we'll use + orientation for all segments
                f.write(f"L\t{src}\t+\t{dst}\t+\t{len(src)-1}M\n")


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
