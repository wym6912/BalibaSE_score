******************************************************************************

	        bali_score evaluation program
                  (version 4.0, Jan 2011)

******************************************************************************

This README contains help with INSTALLATION

Please e-mail bug reports/complaints/suggestions to
	   Julie Thompson at julie@igbmc.fr
 
******************************************************************************

Compilation:
------------

First of all, you need the expat XML parser installed on your machine. If
this is not already done, you can get the parser from http://expat.sourceforge.net/
You should then edit the makefile supplied, changing the EXPAT_INC and EXPAT_LIB paths for your
system.

You make the program with:
make -f makefile

This produces the executable files bali_score bali_score_block bali_score_reliable. 

Usage:
______
To calculate overall column score (CS):
Usage: bali_score ref_aln test_aln 
                where ref_aln       reference alignment in xml format 
                      test_aln      test alignment in msf or fasta format 

To calculate individual block column scores (BCS):
Usage: bali_score_block ref_aln test_aln 
                where ref_aln       reference alignment in xml format 
                      test_aln      test alignment in msf or fasta format 


To calculate overall column score (CS), ignoring sequences with discrepencies (annotated with SEQERR features):
Usage: bali_score_reliable ref_aln test_aln 
                where ref_aln       reference alignment in xml format 
                      test_aln      test alignment in msf or fasta format 


