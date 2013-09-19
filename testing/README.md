# Testing `pareto.py`

These are just manual sanity checks at the moment, based on data I have lying around.

Input data:

* in10.txt:  Ten objectives in columns 27-36 inclusive.  With the epsilons `0.15 30.0 6.0 0.03 30.0 3000.0 150.0 0.3 3.0 0.3`, there should be 632 solutions.
* in3.txt:  Three objectives in columns 27-29 inclusive.  With the epsilons `0.1 0.1 0.1`, there should be 72 solutions.
* in3_{a,b,c,d,e,f}:  Same data as in3.txt, in six input files.  Tests using multiple inputs and getting contribution.

* defective inputs:
    * in3_badcolumn.txt:  Bad data in row 33 (line 34).  Should raise an informative exception when an attempt is made to sort it.
    * in3_shortrow.txt:  Row 4 (line 5) has only 29 columns rather than 30.  Should raise an informative exception when an attempt is made to sort it.
    * in3_longrow.txt  Row 22 (line 23) is double length (as if resulting from an accidental line join in a text editor.)  Right now we ignore the extra columns and pass them through.  Maybe warn about this in future?
* special case inputs:
    * in3_comment.txt  Line 26 has a `#` at the front.  Should produce data identical to three.ref if --comment="#" is specified at the command line, and cause an exception otherwise.
    * in3_blank.txt:  Line 178 is blank.  Should produce data identical to three.ref if --blank is specified at the command line, and cause an exception otherwise.
    * in3_header.txt:  Has a one line header.  Should produce data identical to three.ref if --header=1 is specified at the comand line, and cause an exception otherwise.
    * in3_comment_blank_header.txt:  Has blank lines, comments, and a two-line header.  Should produce data identical to three.ref if appropriate options are specified.
    * in3_obj.txt:  Has only the three objectives columns.  Tests the "assume all columns are objectives" feature.  Also tests the corner case of specifying the wrong number of epsilons when the number of objectives has been determined automatically, if run with a deficient number of epsilons specified.

Output data:

* ten.ref: the 632 solutions that should result from sorting in10.txt.  These were generated using the `ReferenceSetMerger` of `moeaframework`, so they have been externally validated.
* three.ref: the 72 solutions that should result from sorting in3.txt.  These were generated using the `ReferenceSetMerger` of `moeaframework`, so they have been externally validated.

Convenience scripts:

* doit10: reference invocation for sorting in10.txt
* doit3: reference invocation for sorting in3.txt
* comment: reference invocation for sorting in3_comment.txt
* blank: reference invocation for sorting in3_blank.txt
* header: reference invocation for sorting in3_header.txt
* comment_blank_header: reference invocation for sorting in3_comment_blank_header.txt
* contribution: reference invocation for getting contribution information.  Optionally, invoke with --line-number.

Other:

* performance: Totally unscientific performance tracking.  This is just a point of reference for how fast the sort went on a particular machine at a particular time.  It's a baseline to check for major regressions.
