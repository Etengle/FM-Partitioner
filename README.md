# FM_Partitioner
**A 2-way FM partitioner made by Slighten** 
 * A 2-way Partitioning is a technique to divide a circuit (netlist) into 2 parts (called A and B) 
    + The objective is to make the __cut size__ (# of edges) in between A and B __the smaller the better__
    + The constraint is to make the __difference between the size of A and the size of B__ not succeed <br>
      a particular value (e.g. total size divided by 10)
 * FM stands for [Fiduccia-Mattheyses algorithm](https://en.wikipedia.org/wiki/Fiduccia-Mattheyses_algorithm?oldformat=true)
 * For more details you can see [KL (Kernighan–Lin) algorithm](https://en.wikipedia.org/wiki/Kernighan%E2%80%93Lin_algorithm?oldformat=true)
