Welcome to the first version of PCAdapt.
genotype file must contain snps in rows and individuals in column. values are noted with 0, 1, 2 or 9 in case of missing data.

1. To compile type:
	$> make lapack
	$> make

2. To run PCAdapt:
	$> ./PCAdapt
	//will display the manual if compilation went OK
	$> ./PCAdapt -i Example/Data4pops_1 -o My4popsResults -K 5 -b 100 -s 200
	//Will run an example with 400 individuals in 4 pops, and 8000 independent snps, where 400 of them are under selection.

3. To run FastPCAdapt
	$> ./PCAdapt fast
        //will display the manual if compilation went OK
        $> ./PCAdapt fast -i Example/Data4pops_1 Example/Data4pops_2 -o My4popsResults -K 3 
        //Will run an example with 400 individuals in 4 pops, and 10000 independent snps in two different files, where 400 of them are under selection.

4. To clean your folder from executable files:
	$> make clean
	// removes the files in obj and PCAdapt
	$> make realclean
	// removes the files in obj, obj_lapack, and PCAdapt

version 05/26/14:
Computations of Bayes Factors, Posterior Odds, and conditionnal probabilities.
Enhanced SVD initialisation.
Potts model.
Use of different values of Ck.
scaling of the variables by default.
Version of Fast PCAdapt!
