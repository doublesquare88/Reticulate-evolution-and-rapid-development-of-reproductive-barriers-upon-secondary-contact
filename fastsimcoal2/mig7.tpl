//Parameters for the coalescence simulation program : simcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
Npop0
Npop1
Npop2
Npop3
//Haploid samples sizes and samples age 
11
11
11
11
//Growth rates: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
5
//Migration matrix 0
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//Migration matrix 1
0 0 0 0
0 0 0 mig13
0 0 0 0
0 mig13 0 0
//Migration matrix 2
0 0 0 mig03
0 0 0 0
0 0 0 0
mig03 0 0 0
//Migration matrix 3
0 mig01 0 0
mig01 0 0 0
0 0 0 0
0 0 0 0
//Migration matrix 4
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
7 historical event
Tmigstop 0 0 0 1 0 0 
TmigBEA 0 0 0 Resize6 0 1 
TmigBEB 0 0 0 Resize5 0 2 
TmigAM 0 0 0 Resize4 0 3 
Tdiv2 3 2 1 Resize3 0 4
Tdiv2 0 2 1 Resize2 0 4
Tdiv1 2 1 1 Resize1 0 4  
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci