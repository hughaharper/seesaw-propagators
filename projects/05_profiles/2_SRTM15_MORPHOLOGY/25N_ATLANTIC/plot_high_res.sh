#!/usr/bin/env bash

## SOURCE SUBROUTINES
source sub_routines.sh

## SET PLOT VARIABLES
gmt gmtset PS_MEDIA = A4
gmt gmtset PS_PAGE_ORIENTATION = PORTRAIT
title_font="+f12p,Times-Bold+cTC"
pf=high_res_bathy_VGG.ps
ko="-K -O"

## INPUT GRIDS
bathy_grid="../../SRTM15+V2.nc"
VGG_grid="../../Marine_VGG_V27-1.nc"


## SET ATLANTIC REGION
R="-R-49.5/-40.5/23/28"
J="-JY-45/25.5/7c"

## MAKE BATHY CPT
gmt makecpt -CSRTM15+V2.cpt -G-7000/0 > tmp.cpt
gmt makecpt -Ctmp.cpt -T-5000/-500/250 -Z > bathy.cpt


#% START GMT PLOTTING MODE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gmt psxy ${R} ${J} -Y2c -X2c -T -K > ${pf}

## PLOT ATLANTIC BATHYMETRY GRID ---------------------------------------------------------------------------------------

## PLOT BATHY WITH NO OVERLAYS
X=0c
Y=21c
gmt grdcut ${bathy_grid} ${R} -GAtlantic_bathy.nc
@plot_grid Atlantic_bathy.nc bathy.cpt ${R} ${J} ${X} ${Y} I

## PLOT BATHY WITH OVERLAYS
X=0c
Y=-6c
gmt grdcut ${bathy_grid} ${R} -GAtlantic_bathy.nc
@plot_grid Atlantic_bathy.nc bathy.cpt ${R} ${J} ${X} ${Y} I

## REPICK SSP BASED ON DEEPEST POINTS SAMPLED ON 10 km CROSS PROFILES
@auto_pick 012a.xy Atlantic_bathy.nc

## CALCULATE CROSS PROFILES
@generate_profiles 012a.simple Atlantic_bathy.nc
#@generate_profiles 012b.xy Atlantic_bathy.nc

#@generate_profiles 013a.xy Atlantic_bathy.nc
#@generate_profiles 013b.xy Atlantic_bathy.nc

## PLOT EVERY 4th CROSS PROFILE
gmt convert STACKS/012a_profiles.txt -Q0:4:1000 |\
    gmt psxy ${R} ${J} -W0.25p,white ${ko} >> ${pf}

## PLOT CROSS PROFILES
#gmt convert STACKS/012b_profiles.txt -Q0:4:1000 |\
#    gmt psxy ${R} ${J} -W0.25p,white ${ko} >> ${pf}

## PLOT SSPS
gmt psxy 012a.xy ${R} ${J} -W0.5p,white ${ko} >> ${pf}
gmt psxy 012b.xy ${R} ${J} -W0.5p,white ${ko} >> ${pf}


gmt psxy AUTO_PICK/012a_profiles/012a_auto_picked.longlat ${R} ${J} -W0.2p,pink ${ko} >> ${pf}
gmt psxy AUTO_PICK/012a_profiles/012a_auto_picked.filtered ${R} ${J} -W0.2p,black ${ko} >> ${pf}
gmt psxy AUTO_PICK/012a_profiles/012a.simple ${R} ${J} -W0.2p,grey ${ko} >> ${pf}

## PLOT LABELS
font="+f10p,Times-Bold,white"
tail -1 012a.xy | awk '{ print $1, $2, "(c)" }' | gmt pstext ${R} ${J} -F${font} -D-0.25c/0.1c -N ${ko} >> ${pf}
tail -1 012b.xy | awk '{ print $1, $2, "(d)" }' | gmt pstext ${R} ${J} -F${font} -D0.25c/0.1c -N ${ko} >> ${pf}


##----------------------------------------------------------------------------------------------------------------------


## PLOT ATLANTIC VGG GRID ----------------------------------------------------------------------------------------------

## NO OVERLAYS
X=10c
Y=6c
gmt grdcut ${VGG_grid} ${R} -GAtlantic_VGG.nc
@plot_grid Atlantic_VGG.nc Marine_VGG.cpt ${R} ${J} ${X} ${Y}

## WITH OVERLAYS
X=0c
Y=-6c
@plot_grid Atlantic_VGG.nc Marine_VGG.cpt ${R} ${J} ${X} ${Y}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT SCALEBARS
gmt psscale -Dx-9.6c/-1c+w6c/0.2c+jML+h+e -Cbathy.cpt -W0.001 -Bf0.5a1:"Depth (km)": ${ko} >> ${pf}
gmt psscale -Dx0.5c/-1c+w6c/0.2c+jML+h+e -CMarine_VGG.cpt -Bf10a20:"VGG (Etovos)": ${ko} >> ${pf}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT PROFILE 1
Rp="-R-25/25/-0.2/2"
Jp="-JX17c/4c"
gmt psxy ${Rp} ${Jp} -T -X-10C -Y-7c ${ko} >> ${pf}

### PLOT RAW PROFILES
profile_dir="STACKS/012a_profiles/"
g=15
for f in ${profile_dir}/*.xy; do
    gmt psxy ${f} ${Rp} ${Jp} -W0.25p,${g}/${g}/${g} ${ko} >> ${pf}
    ((g++))
    ((g++))
done

## PLOT STAT LINES
gmt psxy ${profile_dir}/relative.std_plot ${Rp} ${Jp} -Gblue@90 ${ko} >> ${pf}
gmt psxy ${profile_dir}/relative.std_plot ${Rp} ${Jp} -W0.5p,blue@60 ${ko} >> ${pf}
gmt psxy ${profile_dir}/relative.mean ${Rp} ${Jp} -W2p,blue ${ko} >> ${pf}
gmt psxy ${profile_dir}/relative.median ${Rp} ${Jp} -W1p,lightblue ${ko} >> ${pf}

## PLOT LABELS
font="+f12p,Times-Bold,black"
echo "-25 2.2 S" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
echo "25 2.2 N" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
echo "-27.5 2.0 (c)" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}

## PLOT BASEMAP
gmt psbasemap ${Rp} ${Jp} -Bf5a10:"Distance (km)":/f0.5a1:"Rel. elevation (km)":nSeW ${ko} >> ${pf}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT PROFILES 2
gmt psxy -T ${Rp} ${Jp} -Y-5.5c ${ko} >> ${pf}

# PLOT RAW PROFILES
profile_dir2="STACKS/012a.simple_profiles"
g=15
for f in ${profile_dir2}/*.xy; do
    gmt psxy ${f} ${Rp} ${Jp} -W0.25p,${g}/${g}/${g} ${ko} >> ${pf}
    ((g++))
    ((g++))
done

## PLOT STAT LINES
echo "PLOTTING STD"
gmt psxy ${profile_dir2}/relative.std_plot ${Rp} ${Jp} -Gred@90 ${ko} >> ${pf}
gmt psxy ${profile_dir2}/relative.std_plot ${Rp} ${Jp} -W0.5p,red@60 ${ko} >> ${pf}
gmt psxy ${profile_dir2}/relative.mean ${Rp} ${Jp} -W2p,red ${ko} >> ${pf}
gmt psxy ${profile_dir2}/relative.median ${Rp} ${Jp} -W1p,lightred ${ko} >> ${pf}

## PLOT LABELS
font="+f12p,Times-Bold,black"
echo "-25 2.2 S" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
echo "25 2.2 N" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}
echo "-27.5 2.0 (d)" | gmt pstext ${Rp} ${Jp} -F${font} -D0.0c/0.0c -N ${ko} >> ${pf}


## PLOT BASEMAP
gmt psbasemap ${Rp} ${Jp} -Bf5a10:"Distance (km)":/f0.5a1:"Rel. elevation (km)":nSeW ${ko} >> ${pf}

##----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#% END GMT
gmt psxy -R -J -O -T >> ${pf}
#% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#% CONVERT TO PNG
gmt psconvert ${pf} -Tf -E150
# osascript -e 'quit app "Preview"' ##% CLOSE PREVIEW
# open ${pf%.ps}.pdf






