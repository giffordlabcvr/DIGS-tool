DIGS screens will often generate redundant, overlapping & fragmented BLAST hits. Fragmented hits can occur because BLAST emphasizes local similarity, and will often fragment a matching region of sequence into several contiguous hits. Of course, this may also occur because a match is genuinely discontiguous.

Overlapping or redundant BLAST hits arise when one or more probes in the query set are closely related, and therefore generate similar - though not necessarily identical - hits. 

The DIGS tool can be configured to deal with these contingencies in different ways, using parameters specified in the control file.

There are three 'redundancy modes' that can be set:

* 1: all BLAST hits are treated as distinct, regardless of whether they overlap at all, are within close range of one another, or are entirely redundant. This mode is only suitable when the aim of using the DIGS tool is simply to establish the presence or absence of specific, relatively rare (i.e. non-repetitive) sequence.

* 2: hits will be merged if they are in the same orientation, have the same value for 'assigned_gene' and are within a specified range of one another.

* 3: hits will be merged if they are in the same orientation, have the same value for 'assigned_gene' and 'assigned_name', and are within a specified range of one another.