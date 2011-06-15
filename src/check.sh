#!/bin/sh
#
#

POLYDIR=$srcdir/src/tests
RESULTSDIR=$srcdir/src/results
UNISOLVE=$srcdir/src/unisolve/unisolve-test

for polynomial_res in $RESULTSDIR/*.res; do
	pol_file=`echo $polynomial_res | sed s/results/tests/ | sed s/\.res/\.pol/`
	echo -n "Solving $pol_file..."
	$UNISOLVE $pol_file $polynomial_res
	if [ "$?" != "0" ]; then
		echo "failed"
		exit 1;
	fi
	echo "done"
done
