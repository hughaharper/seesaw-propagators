#!/usr/bin/env bash

## SET PLOT VARIABLES
gmt gmtset PS_MEDIA = A3
gmt gmtset PS_PAGE_ORIENTATION = PORTRAIT
title_font="+f12p,Times-Bold+cTC"
pf=profiles_in_series.ps
ko="-K -O"


##----------------------------------------------------------------------------------------------------------------------
## PLOT PROFILE 1
Rp="-R-25/25/-0.2/2.2"
Jp="-JX17c/6c"

## START GMT
gmt psxy ${Rp} ${Jp} -T -X2c -Y36c -K > ${pf}

### PLOT RAW PROFILES
ssp=$1
profile_dir=STACKS/${ssp}_profiles/

total=$(echo "110/5" | bc -l)
gmt makecpt -Chaxby -T1/${total}/1 > profiles.cpt

i=1
for p in $(seq 1 5 110); do
#    echo ${i}
    awk -v i=${i} '{ if (NR==1) print "> -Z"i; else print $1, $2}' ${profile_dir}/${ssp}_profile_${p}_relative.xy |\
        gmt psxy ${f} ${Rp} ${Jp} -W2p+cl -Cprofiles.cpt -N ${ko} >> ${pf}
    gmt psxy -T -Y-1.7c ${Rp} ${Jp} ${ko} >> ${pf}
    ((i++))
done

### PLOT STAT LINES
#gmt psxy ${profile_dir}/relative.std_plot ${Rp} ${Jp} -Gblue@90 ${ko} >> ${pf}
#gmt psxy ${profile_dir}/relative.std_plot ${Rp} ${Jp} -W0.5p,blue@60 ${ko} >> ${pf}
#gmt psxy ${profile_dir}/relative.mean ${Rp} ${Jp} -W2p,blue ${ko} >> ${pf}
#gmt psxy ${profile_dir}/relative.median ${Rp} ${Jp} -W1p,lightblue ${ko} >> ${pf}
#
### PLOT LABELS
#font="+f12p,Times-Bold,black"
#echo "-25 2.2 S" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
#echo "25 2.2 N" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
#echo "-27.5 2.0 (c)" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
#
### PLOT BASEMAP
#gmt psbasemap ${Rp} ${Jp} -Bf5a10:"Distance (km)":/f0.5a1:"Rel. elevation (km)":nSeW ${ko} >> ${pf}

#-----------------------------------------------------------------------------------------------------------------------
## END GMT
gmt psxy -R -J -O -T >> ${pf}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## CONVERT TO PNG
gmt psconvert ${pf} -Tf -E150
# osascript -e 'quit app "Preview"' ##% CLOSE PREVIEW
# open ${pf%.ps}.pdf
