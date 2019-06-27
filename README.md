ABC-SeS: ABC framework using random forest coupled with the FDSS or SFS

The main compute_sfs.py script can:
* Compute the FDSS within and between populations 
* Compute the unfolded/folded site frequency spectrum within and between populations
* Compute the number of differences between chr within and between populations
using SNPs simulated with Hudson's ms.

Arguments:

  -h, --help 			show this help message and exit

  -np NPOP, --nrpop NPOP	Number of populations

  -nc NCHR, --nrchr NCHR	Vector containing the number of chr for each pop (es:
                        	2,2,2 for a 3 population comparison)

  -w LLW, --within LLW  	nr of categories that compose the whithin pop freq table 
				(it has to be ajusted based on the expected within population polymorphism)

  -b LLB, --between LLB		nr of categories that compose the between pop freq table
				(it has to be ajusted based on the expected between population polymorphism)

  -s, --segSitesPartition	compute the private and the shared segregating sites between pairs of populations
				instead of pairwise differeces, with a minimum of two pops

  -sfs, --siteFrequencySpectrum	compute the 1D or 2D unfolded sfs instead of pairwise differeces

  -folded, --folded     	fold the site frequency spectrum

  -d, --debug          		print delimiters ("|") between pops and pairwise comparison only using -s




