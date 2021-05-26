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
bathy_grid=SRTM15+V2.nc
VGG_grid=Marine_VGG_V27-1.nc

## SET MH370 REGION
R1="-R89.5/97.5/-37.5/-32.5"
J1="-JY93.5/-35/7c"

## SET ATLANTIC REGION
R2="-R-47/-39/27/32"
J2="-JY-43/29.5/7c"

## MAKE BATHY CPT
gmt makecpt -CSRTM15+V2.cpt -G-7000/0 > tmp.cpt
gmt makecpt -Ctmp.cpt -T-5000/-500/250 -Z > bathy.cpt

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

## CALCULATE CROSS PROFILES
#@generate_profiles 037b.xy mh370_bathy.nc
#@generate_profiles 038b.xy mh370_bathy.nc

@generate_profiles Atlantic_N28.5a.xy Atlantic_bathy.nc
@generate_profiles Atlantic_N28.5b.xy Atlantic_bathy.nc

#@generate_profiles Atlantic_N29a.xy Atlantic_bathy.nc
#@generate_profiles Atlantic_N29b.xy Atlantic_bathy.nc


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#% START GMT PLOTTING MODE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gmt psxy ${R1} ${J1} -Y2c -X2c -T -K > ${pf}

## PLOT MH370 BATHYMETRY GRID ------------------------------------------------------------------------------------------
X=0c
Y=21c
gmt grdcut ${bathy_grid} ${R1} -Gmh370_bathy.nc
@plot_grid mh370_bathy.nc bathy.cpt ${R1} ${J1} ${X} ${Y} I
gmt psxy DZ/037b.xy ${R1} ${J1} -W1p,white ${ko} >> ${pf}
gmt psxy DZ/038b.xy ${R1} ${J1} -W1p,white ${ko} >> ${pf}
##----------------------------------------------------------------------------------------------------------------------


## PLOT MH370 VGG GRID -------------------------------------------------------------------------------------------------
X=10c
Y=0c
gmt grdcut ${VGG_grid} ${R1} -Gmh370_VGG.nc
@plot_grid mh370_VGG.nc Marine_VGG.cpt ${R1} ${J1} ${X} ${Y}
##----------------------------------------------------------------------------------------------------------------------


## PLOT ATLANTIC BATHYMETRY GRID ---------------------------------------------------------------------------------------
X=-10c
Y=-8c
gmt grdcut ${bathy_grid} ${R2} -GAtlantic_bathy.nc
@plot_grid Atlantic_bathy.nc bathy.cpt ${R2} ${J2} ${X} ${Y} I

## PLOT CROSS PROFILES
gmt convert DZ/STACKS/Atlantic_N28.5a_profiles.txt -Q0:4:1000 |\
    gmt psxy ${R2} ${J2} -W0.45p,150/150/150 ${ko} >> ${pf}

## PLOT CROSS PROFILES
gmt convert DZ/STACKS/Atlantic_N28.5b_profiles.txt -Q0:4:1000 |\
    gmt psxy ${R2} ${J2} -W0.45p,150/150/150 ${ko} >> ${pf}

## PLOT SSPS
gmt psxy DZ/Atlantic_N28.5a.xy ${R2} ${J2} -W1p,red ${ko} >> ${pf}
gmt psxy DZ/Atlantic_N28.5b.xy ${R2} ${J2} -W1p,blue ${ko} >> ${pf}

gmt psxy DZ/Atlantic_N29a.xy ${R2} ${J2} -W1p,white ${ko} >> ${pf}
gmt psxy DZ/Atlantic_N29b.xy ${R2} ${J2} -W1p,white ${ko} >> ${pf}


##----------------------------------------------------------------------------------------------------------------------


## PLOT ATLANTIC VGG GRID ---------------------------------------------------------------------------------------
X=10c
Y=0c
gmt grdcut ${VGG_grid} ${R2} -GAtlantic_VGG.nc
@plot_grid Atlantic_VGG.nc Marine_VGG.cpt ${R2} ${J2} ${X} ${Y}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT SCALEBARS
gmt psscale -Dx-9.6c/7c+w6c/0.2c+jML+h+e -Cbathy.cpt -W0.001 -Bf0.5a1:"Depth (km)": ${ko} >> ${pf}
gmt psscale -Dx0.5c/7c+w6c/0.2c+jML+h+e -CMarine_VGG.cpt -Bf10a20:"VGG (Etovos)": ${ko} >> ${pf}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT PROFILE 1
Rp="-R-25/25/-0.2/2"
Jp="-JX17c/4c"
gmt psxy ${Rp} ${Jp} -T -X-10C -Y-6c ${ko} >> ${pf}

## PLOT RAW PROFILES
profile_dir="DZ/STACKS/Atlantic_N28.5b_profiles"
#profile_dir="DZ/STACKS/038b_profiles"
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

## PLOT BASEMAP
gmt psbasemap ${Rp} ${Jp} -Bf5a10:"Distance (km)":/f0.5a1:"Rel. elevation (km)":nSeW ${ko} >> ${pf}
##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## PLOT PROFILES 2
gmt psxy -T ${Rp} ${Jp} -Y-5.5c ${ko} >> ${pf}

## PLOT RAW PROFILES
profile_dir2="DZ/STACKS/Atlantic_N28.5a_profiles"
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







