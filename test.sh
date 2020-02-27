#! /bin/bash

source parameterfile_Ref.init

if grep -q "chr" $OUTDIR/$PROJECT/$PROJECT.input
then

echo "running script for chr present"
perl 01_DelP_findcorrespondinginsertion_v3.3.pl -t $RM_TRACK -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/RM_intervals.out -te $1

else

echo "script for chr absent"
perl 01_DelP_findcorrespondinginsertion_v3.3.pl -t <(sed 's/chr//g' <(cat $RM_TRACK)) -f $OUTDIR/$PROJECT/$PROJECT.input -p $OUTDIR/$PROJECT/RM_intervals.out -te $1

fi
