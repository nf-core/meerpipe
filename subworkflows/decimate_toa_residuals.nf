

// Convert nchan and npols to lists
nchans = params.nchans.split(',').collect { it.toInteger() }
npols  = params.npols.split(',').collect  { it.toInteger() }


process GRAB_ALL_PAIRS {
    label 'ephem_template'

    input:
    val pulsar

    output:
    path "all_pairs.csv"

    """
    grab_all_pairs ${pulsar}
    """
}


process DECIMATE {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/decimated", mode: 'copy', pattern: "${pulsar}_${utc}_zap.*.ar"

    input:
    tuple val(pulsar), val(utc), val(beam), val(dur), val(pipe_id), path(cleaned_archive), val(snr)

    output:
    tuple val(pulsar), val(utc), val(beam), val(dur), val(pipe_id), path("${pulsar}_${utc}_zap.*.ar")

    """
    for nchan in ${nchans.join(' ')}; do
        if ${params.use_max_nsub}; then
            # Calculate nsub to get desired TOA S/N
            max_nsub=\$(python -c "import math; print(math.floor(1/\$nchan * (${snr}/${params.tos_sn}) ** 2))")

            input_nsub=\$(vap -c nsub ${cleaned_archive} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
            if [ \$max_nsub -gt \$input_nsub ]; then
                # Greater than input nsub so set input as max
                max_nsub=\$input_nsub
            fi
            if [ \$max_nsub -eq 0 ]; then
                # Not enough SN so only make a fully time scrunched
                nsubs="1"
            else
                nsubs="1 \$max_nsub"
            fi
        else
            nsubs="1"
        fi

        # Make a max_nsub decimation and a time scrunched decimation
        for nsub in \$nsubs; do
            # Make full stokes and/or polarisation scrunched
            for stokes in ${npols.join(' ')}; do
                if [ \${stokes} -eq 1 ]; then
                    # Polarisation scrunch option
                    stokes_op="-p"
                else
                    stokes_op=""
                fi

                echo "Decimate nsub=\${nsub}  nchan=\${nchan} stokes=\${stokes}"
                pam --setnsub \${nsub} --setnchn \${nchan} -S \${stokes_op} -e \${nchan}ch\${stokes}p\${nsub}t.ar ${cleaned_archive}
            done
        done
    done
    """
}



process GENERATE_TOAS {
    label 'cpu'
    label 'psrchive'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/timing/${project_short}", mode: 'copy', pattern: "*.{residual,tim,par,std}"

    input:
    tuple val(pulsar), val(project_short), path(ephemeris), path(template), val(utc), val(beam), val(dur), val(pipe_id), path(decimated_archives)

    output:
    tuple val(pulsar), val(project_short), path(ephemeris), path(template), val(pipe_id), path("*.tim"), path("*.residual")

    """
    # Loop over each DECIMATEd archive
    for ar in ${decimated_archives.join(' ')}; do
        if [[ \$ar == *"ch4p"* ]]; then
            # Skip if it is a full Stokes archive
            continue
        fi

        # Grab archive nchan and nsub
        nchan=\$(vap -c nchan \$ar | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        nsub=\$( vap -c nsub  \$ar | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        # Grab template nchan
        tnchan=\$(vap -c nchan ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)

        # Use portrait mode if template has more frequency channels
        if [ "\$tnchan" -gt "\$nchan" ]; then
            port="-P"
        else
            port=""
        fi

        echo "Generating TOAs for \${ar}.tim\n----------------------------------"
        pat -jp \$port  -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM \$ar  > \${ar}.tim
    done

    # Create residuals for time largest archive
    largest_archive=\$(ls ${pulsar}_${utc}_zap.${nchans.max()}ch1p*t.ar | tail -n 1)
    tempo2_wrapper.sh \${largest_archive} ${ephemeris}
    """
}
    // DM corrected stuff
    //     echo "Correct for DM\n----------------------------------"
    //     dm=\$(grep DM ${dm_results} | cut -d ' ' -f 2)
    //     pam -D \$dm -e ar.dm_corrected \$ar

    //     echo "Generating TOAs for DM corrected archive\n----------------------------------"
    //     pat -jp \$port  -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM \$ar.dm_corrected  > \${ar}.dm_corrected.tim
    // # And largest DM corrected archive
    // tempo2_wrapper.sh \${largest_archive}.dm_corrected ${ephemeris}



process UPLOAD_TOAS {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(pulsar), val(project_short), path(ephemeris), path(template), val(pipe_id), path(toas), path(residuals)

    output:
    tuple val(pulsar), val(project_short), path(ephemeris)


    """
    #!/usr/bin/env python

    import os
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_graphql_id, get_rest_api_id
    from psrdb.tables.toa import Toa
    from psrdb.tables.template import Template

    import psrdb
    print(psrdb.__file__)

    logger = setup_logging(console=True, level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    toa_client      = Toa(client)
    template_client = Template(client)

    # Upload template (if not uploaded already)
    template = os.path.realpath("${template}")
    template_band = template.split("/")[-3]
    template_project = template.split("/")[-2]
    template_response = template_client.create(
        "${pulsar}",
        template_band,
        template,
        project_short=template_project,
    )
    logger.debug(template_response)
    template_id = get_rest_api_id(template_response, logging.getLogger(__name__))

    # Upload TOAs
    for toa_file in ["${toas.join('","')}"]:
        if "dm_corrected" in toa_file:
            dmcorrected = True
            # TODO don't skip this if it's useful
            continue
        else:
            dmcorrected = False

        # Work out if this is the minimum or maximum number of subints
        nchan = toa_file.split("_zap.")[-1].split("ch")[0]
        nsub = toa_file.split("_zap."+nchan+"ch1p")[-1].split("t.ar")[0]
        toas_same_nchan = glob("*_zap." + nchan + "ch*" + toa_file.split(".ar.")[1])
        nsubs_list = []
        for toa in toas_same_nchan:
            nsubs_list.append(int(toa.split("_zap."+nchan+"ch1p")[-1].split("t.ar")[0]))
        minimum_nsubs = False
        maximum_nsubs = False
        if max(nsubs_list) == int(nsub):
            maximum_nsubs = True
        if min(nsubs_list) == int(nsub):
            minimum_nsubs = True

        logger.info(f"Uploading Toa file {toa_file} with maximum_nsubs={maximum_nsubs} and minimum_nsubs={minimum_nsubs}")
        with open(toa_file, "r") as f:
            toa_lines = f.readlines()
            toa_response = toa_client.create(
                ${pipe_id},
                "${project_short}",
                "${ephemeris}",
                template_id,
                toa_lines,
                dmcorrected,
                minimum_nsubs,
                maximum_nsubs,
            )
            if toa_response.status_code not in (200, 201):
                logger.error("Failed to upload TOA")
                exit(1)
            logger.info(get_graphql_id(toa_response, "toa", logger))
    """
}



process GENERATE_RESIDUALS {
    label 'psrdb_tempo2'

    maxForks 1

    input:
    tuple val(pulsar), val(project_short), path(ephemeris)

    """
    # Loop over each of the TOA filters
    for min_or_max_sub in "--minimum_nsubs" "--maximum_nsubs"; do
        for nchan in ${nchans.join(' ')}; do
            echo -e "\\nDownload the toa file and fit the residuals for \${min_or_max_sub#--} \${nchan}\\n--------------------------\\n"
            psrdb toa download ${pulsar} --project ${project_short} \$min_or_max_sub --nchan \$nchan
            echo -e "\\nGenerating residuals for \${min_or_max_sub#--} \${nchan}\\n--------------------------\\n"
            tempo2_wrapper.sh toa_${pulsar}_\${min_or_max_sub#--}_nchan\${nchan}.tim ${ephemeris}
            if [ -f "toa_${pulsar}_\${min_or_max_sub#--}_nchan\${nchan}.tim.residual" ]; then
                echo -e "\\nUpload the residuals\\n--------------------------\\n"
                psrdb residual create toa_${pulsar}_\${min_or_max_sub#--}_nchan\${nchan}.tim.residual
            fi
        done
    done
    """
}

workflow DECIMATE_TOA_RESIDUALS {
    take:
        files_and_meta // channel

    main:
        // Grab all ephemeris and template pairs for each pulsar
        GRAB_ALL_PAIRS( files_and_meta.groupTuple().map {
                pulsar, utc, beam, dur, pipe_id, cleaned_archive, snr ->
                pulsar
            }
        )

        DECIMATE( files_and_meta )

        // Decimate into different time and freq chunks using pam
        // DECIMATE( files_and_meta )

        // // Generate TOAs
        GENERATE_TOAS(
            GRAB_ALL_PAIRS.out.splitCsv().combine(DECIMATE.out, by: 0)
        )

        if ( params.upload ) {
            UPLOAD_TOAS( GENERATE_TOAS.out )

            // For each pulsar (not each obs), download all toas and fit residuals
            GENERATE_RESIDUALS(
                UPLOAD_TOAS.out.groupTuple().cross( GRAB_ALL_PAIRS.out.splitCsv() ).map{ it[1] }.map {
                    pulsar, project_short, ephem, template -> [ pulsar, project_short, ephem ]
                }
            )
        }
}