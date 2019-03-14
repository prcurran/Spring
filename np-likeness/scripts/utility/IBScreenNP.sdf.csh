#!/bin/csh -f
# usual.sdf.csh from http://zinc.docking.org
# Run this script to download ZINC
# Requires curl (default) or wget, http://wget.docking.org
#
# Thus, to run this script
#         using curl, do:     csh usual.sdf.csh
#         using wget, do:     csh usual.sdf.csh wget
#
setenv base http://zinc.docking.org/db/byvendor/ibsnp
setenv fn .zinc.$$
cat <<+ > $fn
ibsnp_p0.0.sdf.gz
ibsnp_p1.0.sdf.gz
+
if ($#argv>0) then
     wget --base=$base -i < $fn
else
     foreach i (`cat $fn`)
          curl --url $base/$i -o $i
     end
endif
rm -f $fn
# File created on  Sun Jun 21 11:02:30 PDT 2015
# This is the end of the csh script.
