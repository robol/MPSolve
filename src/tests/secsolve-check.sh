
#!/bin/sh
#
#

POLYDIR=$srcdir/src/tests/secsolve
RESULTSDIR=$srcdir/src/results/secsolve
SECSOLVE=$srcdir/src/secsolve/secsolve-test
TMPFILE=$srcdir/secsolve_output

top_dir=`pwd`
cd $RESULTSDIR

passed=0

for precision in 15 50; do
    for polynomial_res in *.res; do
            pol_file=`echo $polynomial_res | sed s/\.res/\.pol/`
            printf "Testing resolution of %15s... " $polynomial_res
            $top_dir/$SECSOLVE $top_dir/$POLYDIR/$pol_file \
              $top_dir/$RESULTSDIR/$polynomial_res -tf -o$precision 2> $top_dir/$TMPFILE
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

    for polynomial_res in *.res; do
            pol_file=`echo $polynomial_res | sed s/\.res/\.pol/`
            printf "Testing resolution of %15s... " $polynomial_res
            $top_dir/$SECSOLVE $top_dir/$POLYDIR/$pol_file \
              $top_dir/$RESULTSDIR/$polynomial_res -td -o$precision 2> $top_dir/$TMPFILE
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
done

exit $passed
