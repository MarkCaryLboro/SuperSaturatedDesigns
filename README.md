# SuperSaturatedDesigns

Super saturated design code

1. Factor screeening for large numbers of factors.
2. Main effects model only.
3. More factors than runs.

Implements optimal supersaturated design algorithms. Supported algorithms are Bayesian (Jones et al), Li and Marley.
The Bayesian method is much faster and as robust as the other two. Consequently, this is highly recommended.

All the algorithms make use of a columnwise pairwise (CWPW) exchange optimal search, which does not require the definition
of a candidate set. With large numbers of factors this is advantageous, as the candidate set will require significant storage.
Also, it is inefficient with a Federov style row-exchange algorithm to search the entire candidate list to implement
the optimal coordinate exchange. Column exchanges require far less computational effort as the number of possible exchanges is 
greatly reduced.

The user notes are contained in file "Bayesian Supersaturated Designs for Factor Screening2_2.docx" and contain background 
information, installation and user instructions. An example is also worked through to demonstrate code use.

The powerpoint presentation describes the Bayesian algorithm and sets it in the context of screening lithium-ion battery variables.

Mark Cary 


 