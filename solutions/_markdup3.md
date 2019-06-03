Picard and samtools uses the alignment positions:
 
   - Both 5' ends of both reads need to have the same positions. 
   - Each reads have to be on the same strand as well.

Another method is to use a kmer approach:
  
   - take a part of both ends of the fragment
   - build a hash table 
   - count the similar hits

Brute force, compare all the sequences.
