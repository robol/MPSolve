
#!/bin/bash
#
#
# INITIAL CONFIGURATION
POLYDIR=$srcdir/src/tests/unisolve
RESULTSDIR=$srcdir/src/results/unisolve
UNISOLVE=$srcdir/src/unisolve/unisolve-test
TMPFILE=$srcdir/unisolve_output

top_dir=`pwd`
cd $RESULTSDIR

passed=0

for precision in 50 100; do
    echo "Checking polynomials asking MPSolve to compute"
    echo "$precision exact digits."
    echo ""

    for polynomial_res in *.res; do
	pol_file=`echo $polynomial_res | sed s/\.res/\.pol/`
	printf "Testing resolution of %15s... " $polynomial_res
	$top_dir/$UNISOLVE $top_dir/$POLYDIR/$pol_file \
	    $top_dir/$RESULTSDIR/$polynomial_res -o$precision 2> $top_dir/$TMPFILE
	if [ "$?" != "0" ]; then
	    echo -e "\033[31;1mfailed\033[0m"
	    echo -e "\n\033[1mError reported:\033[0m"
	    cat < $top_dir/$TMPFILE
	    echo ""
	    passed=1;
	else
		echo -e "\033[32;1mpassed\033[0m"
	fi
	rm -f $top_dir/$TMPFILE
    done

    echo ""
    echo "Test finished"
    echo ""
done

exit $passed
