
#!/bin/sh
#
#

POLYDIR=$srcdir/src/tests/secsolve
RESULTSDIR=$srcdir/src/results/secsolve
SECSOLVE=$srcdir/src/secsolve/secsolve-test
TMPFILE=$srcdir/secsolve_output

top_dir=`pwd`
cd $RESULTSDIR
for polynomial_res in *.res; do
	pol_file=`echo $polynomial_res | sed s/\.res/\.pol/`
	printf "Testing resolution of %15s... " $polynomial_res
	$top_dir/$SECSOLVE $top_dir/$POLYDIR/$pol_file \
	  $top_dir/$RESULTSDIR/$polynomial_res f 2> $top_dir/$TMPFILE
	if [ "$?" != "0" ]; then
		echo -e "\033[31;1mfailed\033[0m"
		echo -e "\n\033[1mError reported:\033[0m"
		cat < $top_dir/$TMPFILE
		echo ""
		exit 1;
	fi
	echo -e "\033[32;1mpassed\033[0m"
	rm -f $top_dir/$TMPFILE
done

for polynomial_res in *.res; do
	pol_file=`echo $polynomial_res | sed s/\.res/\.pol/`
	printf "Testing resolution of %15s... " $polynomial_res
	$top_dir/$SECSOLVE $top_dir/$POLYDIR/$pol_file \
	  $top_dir/$RESULTSDIR/$polynomial_res d 2> $top_dir/$TMPFILE
	if [ "$?" != "0" ]; then
		echo -e "\033[31;1mfailed\033[0m"
		echo -e "\n\033[1mError reported:\033[0m"
		cat < $top_dir/$TMPFILE
		echo ""
		exit 1;
	fi
	echo -e "\033[32;1mpassed\033[0m"
	rm -f $top_dir/$TMPFILE
done

