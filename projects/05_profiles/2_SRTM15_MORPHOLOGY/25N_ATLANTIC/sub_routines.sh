#!/usr/bin/env bash

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## GMT GRID PLOTTING SUBROUTINE
@plot_grid () {
    # SET INPUTS
    grid=$1
    cpt=$2
    R=$3
    J=$4
    X=$5
    Y=$6

    ## PLOT GRID
    if [ ! -z "$7" ]; then
        ## USE GRADIENT
        gmt grdimage ${grid} ${R} ${J} -I+A315+nt0.9 -C${cpt} -X${X} -Y${Y} \
            ${ko} >> ${pf}
    else
        ## DON'T USE GRADIENT
        gmt grdimage ${grid} ${R} ${J} -C${cpt} -X${X} -Y${Y} ${ko} >> ${pf}
    fi
    ## PLOT COAST
    gmt pscoast ${R} ${J} -W0.25p -A5000 -Di ${ko} >> ${pf}

    ## PLOT BASEMAP
    gmt psbasemap ${R} ${J} -Bf1a2:"":/f1a2:"":nSeW ${ko} >> ${pf}
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## CROSS PROFILE SUBROUTINE
@generate_profiles () {
    ## SET INPUTS
    ssp=$1
    grid=$2

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 1.0 GENERATE CROSS PROFILES USING GRDTRACK (-C IS KEY PARAMETER)
    gmt grdtrack ${ssp} -G${grid} -C50k/0.5k/2k+v \
        -Sm+sSTACKS/${ssp%.xy}_stack.txt > STACKS/${ssp%.xy}_profiles.txt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 2.0 EXTRACT INDIVIDUAL PROFILES FROM GRDTRACK OUTPUT
    ssp=${ssp%.xy}
    mkdir STACKS/${ssp}_profiles
    i=0
    awk -v i=${i} -v ssp=${ssp} '{ if ($1 == ">")
        close ("STACKS/"ssp"_profiles/"ssp"_profile_"i".txt") i++;
            else
        print $0 > ("STACKS/"ssp"_profiles/"ssp"_profile_"i".txt")
    }' STACKS/${ssp%.xy}_profiles.txt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 3.0 PROCESS INDIVIDUAL PROFILES (CALCULATE DEPTHS RELATIVE TO DEEPEST
    ##                                  POINT ALONG THE PROFILE)

    ## SET WORKING DIR
    wrk_dir="STACKS/${ssp}_profiles"

    ## EXTRACT X PROFILE VALUES (USED LATER IN SCRIPT)
    awk '{ print $3 }' ${wrk_dir}/*_1.txt > ${wrk_dir}/x.all

    ## RUN PROCESSING LOOP
    i=1
    for f in ${wrk_dir}/*.txt; do
        ## 3.1 FIND DEEPEST POINT IN PROFILE
        deepest_point=$(sort -nk 5 ${f} | head -n 1 | awk '{ print $5}')

        ## 3.2 CALCULATE DEPTHS RELATIVE TO DEEPEST POINT
        awk -v dp=${deepest_point} '{
            print $3, ($5-dp)*0.001 }' ${f} > ${f%.txt}_relative.xy

        ## 3.3 PRINT ALL Y VALUES TO relative.all FOR MEDIAN CALCULATION
        awk -v dp=${deepest_point} '{ print ($5-dp)*0.001 }' ${f} > relative_y

        if [[ ${i} -eq 1 ]]; then
            cat relative_y > ${wrk_dir}/relative.all
        else
            paste ${wrk_dir}/relative.all relative_y > tmp
            mv tmp ${wrk_dir}/relative.all
        fi

        ((i++))
    done
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 4.0 CREATE MEDIAN PROFILE

    ## REMOVE OLD FILE
    rm ${wrk_dir}/median.all

    ## GET LINE COUNT
    line_count=$(wc -l ${wrk_dir}/relative.all | awk '{print $1}')

    ## RUN LOOP
    for lc in $(seq 1 1 ${line_count}); do

        ## EXTRACT ALL VALUES FROM ROW lc AND SORT NUMERICALLY
        awk -v lc=${lc} 'NR==lc''{
                for (i=1;i<=NF;i++) if (i<=NF) print $i
                    }' ${wrk_dir}/relative.all | sort -k1n > tmp

        ## GET MIDDLE LINE
        median_line_number=$(wc -l tmp | awk '{printf "%i", $1/2}')

        ## GET VALUE OF MIDDLE LINE
        median=$(awk -v m=${median_line_number} 'NR==m''{print $1}' tmp)

        ## APPEND MEDIAN VALUE TO MASTER FILE
        echo ${median} >> ${wrk_dir}/median.all
    done
    paste ${wrk_dir}/x.all ${wrk_dir}/median.all > ${wrk_dir}/relative.median
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 5.0 CREATE MEAN PROFILE
    awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i};
            sum/=NF; print sum  }' ${wrk_dir}/relative.all > ${wrk_dir}/mean.all

    paste ${wrk_dir}/x.all ${wrk_dir}/mean.all > ${wrk_dir}/relative.mean
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 6.0 CALCULATE STD
    rm ${wrk_dir}/relative.std

    for lc in $(seq 1 1 ${line_count}); do

        ## GET MEAN VALUE
        mean=$(awk -v lc=${lc} 'NR==lc''{print $2}' ${wrk_dir}/relative.mean)

        ## CALCULATE: [(VALUE-mean)^2]  FOR ALL VALUES OF POINT lc AND THE SUM
        sum=$(awk -v lc=${lc} -v mean=${mean} 'NR==lc''{
            for (i=1;i<=NF;i++) if (i<=NF) print ($i-mean)^2
                }' ${wrk_dir}/relative.all | awk '{sum+=$1} END{print sum}')

        ## GET NUMBER OF PROFILES-1
        N=$(awk 'NR==1''{ print (NF-1)}' ${wrk_dir}/relative.all)

        ## CALCULATE STD DEV FOR THE PROFILE POINT lc
        std=$(echo "sqrt(${sum}/${N})" | bc -l)

        if [[ ${lc} -eq 1 ]]; then
            awk -v lc=${lc} -v std=${std} 'NR==lc''{
                print $1, $2+std, $2-std
                    }' ${wrk_dir}/relative.mean > ${wrk_dir}/relative.std
        else
            awk -v lc=${lc} -v std=${std} 'NR==lc''{
                print $1, $2+std, $2-std
                    }' ${wrk_dir}/relative.mean >> ${wrk_dir}/relative.std
        fi

    done

    awk '{ print $1, $2 }' ${wrk_dir}/relative.std > ${wrk_dir}/relative.std_plot
    awk '{ print $1, $3 }' ${wrk_dir}/relative.std | tail -r >> ${wrk_dir}/relative.std_plot
    awk 'NR==1''{ print $1, $2 }'  ${wrk_dir}/relative.std >> ${wrk_dir}/relative.std_plot

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## CROSS PROFILE SUBROUTINE
@auto_pick () {
    ## SET INPUTS
    ssp=$1
    grid=$2

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 1.0 GENERATE CROSS PROFILES USING GRDTRACK (-C IS KEY PARAMETER)
    gmt grdtrack ${ssp} -G${grid} -C10k/0.5k/2k+v \
        -Sm+sAUTO_PICK/${ssp%.xy}_stack.txt > AUTO_PICK/${ssp%.xy}_profiles.txt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 2.0 EXTRACT INDIVIDUAL PROFILES FROM GRDTRACK OUTPUT
    ssp=${ssp%.xy}
    mkdir AUTO_PICK/${ssp}_profiles
    i=0
    awk -v i=${i} -v ssp=${ssp} '{ if ($1 == ">")
        close ("AUTO_PICK/"ssp"_profiles/"ssp"_profile_"i".txt") i++;
            else
        print $0 > ("AUTO_PICK/"ssp"_profiles/"ssp"_profile_"i".txt")
    }' AUTO_PICK/${ssp%.xy}_profiles.txt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 3.0 PROCESS INDIVIDUAL PROFILES (CALCULATE DEPTHS RELATIVE TO DEEPEST
    ##                                  POINT ALONG THE PROFILE)

    ## SET WORKING DIR
    wrk_dir="AUTO_PICK/${ssp}_profiles"

    ## RUN PROCESSING LOOP
    i=1
    > ${wrk_dir}/${ssp}_auto_picked.longlat
    for f in ${wrk_dir}/*.txt; do
        ## 3.1 FIND DEEPEST POINT IN PROFILE
        deepest_point=$(sort -nk 5 ${f} | head -n 1 | awk '{ print $5}')

        ## 3.2 REPICK THE SSP BASED ON DEEPEST POINTS
        sort -nk 5 ${f} | head -n 1 | awk '{ print $1, $2}' >> ${wrk_dir}/${ssp}_auto_picked.tmp

        ((i++))
    done

    sort -nk 1 ${wrk_dir}/${ssp}_auto_picked.tmp > ${wrk_dir}/${ssp}_auto_picked.longlat

    gmt sample1d ${wrk_dir}/${ssp}_auto_picked.longlat -Fa -T0.5k -fg |\
        filter1d -FC0.1 -E > ${wrk_dir}/${ssp}_auto_picked.filtered

    gmt simplify ${wrk_dir}/${ssp}_auto_picked.filtered -T6k > ${wrk_dir}/${ssp}.simple

    rm ${wrk_dir}/${ssp}_auto_picked.tmp
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------