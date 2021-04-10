#!/bin/bash
# Run this file to generate the figures
# The contour plots take a long time so be patient.

make iondenofphi
make contribs
./contribs
mv plot0001.ps contribs.ps
epstopdf contribs.ps

./iondenofphi -E -c
mv plot0001.ps twobeamf.ps
epstopdf twobeamf.ps
mv plot0002.ps dFdxvb.ps
epstopdf dFdxvb.ps
mv plot0003.ps phinofx.ps
epstopdf phinofx.ps
mv plot0005.ps phinofx2.ps
epstopdf phinofx2.ps
mv plot0006.ps phinofx3.ps
epstopdf phinofx3.ps

#./iondenofphi
./iondenofphi -E -M -U -p.02 -T1 -t.2 -d.2 -c
mv plot0002.ps vsthreshplot.ps
epstopdf vsthreshplot.ps
./iondenofphi -E -M -U -p.5 -T1 -t.2 -d.2 -c
mv plot0002.ps vsthreshplot2.ps
epstopdf vsthreshplot2.ps

# This file requires other code. Don't delete it.
#./chiofv -v.8,1.2,1.,1. -v.2,-1.5,.3,.3 -T4 -x.7 -d
#mv plot0001.ps ionstabplot.ps
#ps2png '-r300 ionstabplot.ps' ionstabplot.pdf

# ion-ion T-threshold contours 20x20 only
./iondenofphi -M -p.5 -c -C 2>/dev/null
mv plot0401.ps contpsi5.ps
epstopdf  contpsi5.ps

./iondenofphi -M -p.1 -c -C 2>/dev/null
mv plot0401.ps contpsi1.ps
epstopdf  contpsi1.ps

./iondenofphi -M -p.5 -i2 -m3 -c
mv plot0001.ps shapes.ps
epstopdf shapes.ps

rm *.ps
