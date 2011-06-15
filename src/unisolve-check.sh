#!/bin/sh
#
#

POLYDIR=$srcdir/src/tests
RESULTSDIR=$srcdir/src/results
UNISOLVE=$srcdir/src/unisolve/unisolve-test

for polynomial_res in $RESULTSDIR/*.res; do
	pol_file=`echo $polynomial_res | sed s/results/tests/ | sed s/\.res/\.pol/`
	printf "Testing resolution of %35s...\t" $pol_file
	$UNISOLVE $pol_file $polynomial_res
	if [ "$?" != "0" ]; then
		echo "\033[34;1mfailed\033[0m"
		exit 1;
	fi
	echo "\033[32;1mpassed\033[0m"
done
