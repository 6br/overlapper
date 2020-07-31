# overlapper

Extract the number of bridges between contigs from long-read alignments to visualize the connection (under development).

### Usage

The input is PAF alignments of long reads against contigs(assembly).

```
stack build
stack exec overlapper-exe --RTS -- /*THE NAME OF PAF*/ /* output json */ +RTS -N12 > /*log file*/
```

The output is the pair of contigs and read counts. 
