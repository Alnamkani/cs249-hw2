To run the overlap graph generator:

`python olc.py -f sequences.fastq -n 4`

To run De Bruijn graph:

`python dbg.py -f sequences.fastq -k 4`

To export the De Bruijn graph as a GFA file:

`python dbg.py -f sequences.fastq -k 4 -g graph.gfa`

To change the output file name:

For OLC: `python olc.py -f sequences.fastq -n 4 -o custom_output.fasta`

For DBG: `python dbg.py -f sequences.fastq -k 4 -o custom_output.fasta`