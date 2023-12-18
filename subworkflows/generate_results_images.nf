

process DM_RM_CALC {
    label 'cpu'
    label 'meerpipe'

    input:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr)

    output:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path("${pulsar}_${utc}_dm_rm_fit.txt")

    script:
    if ( task.attempt > 2 )
        """
        echo "DM: None"     >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "ERR: None"    >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "EPOCH: None"  >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "CHI2R: None"  >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "TRES: None"   >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM: None"     >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM_ERR: None" >> ${pulsar}_${utc}_dm_rm_fit.txt
        """
    else if ( Float.valueOf(snr) > 20.0 )
        """
        echo -e "\\nCreate a max channel archive\\n----------------------------------"
        # Calculate nchan to get desired TOA S/N and make sure it is a factor of the archive channels
        arnchan=\$(vap -c nchan ${cleaned_archive} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        nchan=\$(python -c "import math; raw_nchan=math.floor( (${snr}/10.) ** 2); print(next((factor for factor in range(\$arnchan + 1, 2, -1) if \$arnchan % factor == 0 and factor < raw_nchan and factor <= 64)))")
        if [ \$nchan -gt 16 ]; then
            nchan=16
        fi
        pam --setnchn \${nchan} -T -S -p -e dmcalc ${cleaned_archive}
        pam --setnchn \${nchan} -T -S    -e rmcalc ${cleaned_archive}

        echo -e "\\nCreate TOAs with max channel archive\\n----------------------------------"
        # Grab template nchan
        tnchan=\$(vap -c nchan ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        # Use portrait mode if template has more frequency channels
        if [ "\$tnchan" -gt "\$nchan" ]; then
            port="-P"
        else
            port=""
        fi
        pat -jp \$port -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM ${pulsar}_${utc}_zap.dmcalc > dm.tim

        echo -e "\\nCalc DM with tempo2\\n----------------------------------"
        # Remove dm derivatives
        sed '/^DM[1-9]/d' ${ephemeris} > ${ephemeris}.dm
        echo "MODE 1" >>  ${ephemeris}.dm
        # Fit for DM
        tempo2 -nofit -fit DM -set START 40000 -set FINISH 99999 -f ${ephemeris}.dm -outpar ${ephemeris}.dmfit dm.tim

        input_rm=\$(vap -c rm ${pulsar}_${utc}_zap.rmcalc | tail -n 1| tr -s ' ' | cut -d ' ' -f 2)
        lower_rm=\$(echo "\$input_rm - 34" | bc -l)
        higher_rm=\$(echo "\$input_rm + 34" | bc -l)
        echo -e "\\nFit for RM from \$lower_rm - \$higher_rm \\n----------------------------------"
        rmfit -D -m \$lower_rm,\$higher_rm,200 ${pulsar}_${utc}_zap.rmcalc -K /PNG > rmfit_output.txt

        echo -e "\\nGrab the outputs and write it to a file\\n----------------------------------"
        DM=\$(grep    "^DM "      ${ephemeris}.dmfit | awk '{print \$2}')
        ERR=\$(grep   "^DM "      ${ephemeris}.dmfit | awk '{print \$4}')
        EPOCH=\$(grep "^DMEPOCH " ${ephemeris}.dmfit | awk '{print \$2}')
        CHI2R=\$(grep "^CHI2R "   ${ephemeris}.dmfit | awk '{print \$2}')
        TRES=\$(grep  "^TRES "    ${ephemeris}.dmfit | awk '{print \$2}')
        rm_results=\$(grep "Best RM is" rmfit_output.txt | cut -d ':' -f 2)
        RM=\$(echo     \$rm_results | cut -d '/' -f 1 | cut -d ' ' -f 1)
        RM_ERR=\$(echo \$rm_results | cut -d '/' -f 2 | cut -d ' ' -f 2)

        echo "DM: \${DM}"         >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "ERR: \${ERR}"       >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "EPOCH: \${EPOCH}"   >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "CHI2R: \${CHI2R}"   >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "TRES: \${TRES}"     >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM: \${RM}"         >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM_ERR: \${RM_ERR}" >> ${pulsar}_${utc}_dm_rm_fit.txt
        """
    else
        """
        pam -T -S -p -e dmcalc ${cleaned_archive}
        pdmp -g ${cleaned_archive}.ps/cps -mc 16 *dmcalc

        # Grab the outputs and write it to a file
        DM=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 5)
        ERR=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 6)
        EPOCH=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 2)
        CHI2R=None
        TRES=None

        echo "DM: \${DM}"       >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "ERR: \${ERR}"     >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "EPOCH: \${EPOCH}" >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "CHI2R: \${CHI2R}" >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "TRES: \${TRES}"   >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM: None"         >> ${pulsar}_${utc}_dm_rm_fit.txt
        echo "RM_ERR: None"     >> ${pulsar}_${utc}_dm_rm_fit.txt
        """
}

process GENERATE_IMAGE_RESULTS {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/images", mode: 'copy', pattern: "{c,t,r}*png"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/scintillation", mode: 'copy', pattern: "*dynspec*"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: 'copy', pattern: "results.json"

    input:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path(dm_results)

    output:
    tuple val(pulsar), val(pipe_id), path("*.png"), path("*.dat"), path("*dynspec"), path("results.json")


    """
    # psrplot images
    for i in "raw ${raw_archive}" "cleaned ${cleaned_archive}"; do
        set -- \$i
        type=\$1
        file=\$2
        # Do the plots for raw file then cleaned file
        psrplot -p flux -jFTDp -jC                          -g 1024x768 -c above:l= -c above:c="Stokes I Profile (\${type})"     -D \${type}_profile_fts.png/png \$file
        psrplot -p Scyl -jFTD  -jC                          -g 1024x768 -c above:l= -c above:c="Polarisation Profile (\${type})" -D \${type}_profile_ftp.png/png \$file
        psrplot -p freq -jTDp  -jC                          -g 1024x768 -c above:l= -c above:c="Phase vs. Frequency (\${type})"  -D \${type}_phase_freq.png/png  \$file
        psrplot -p time -jFDp  -jC                          -g 1024x768 -c above:l= -c above:c="Phase vs. Time (\${type})"       -D \${type}_phase_time.png/png  \$file
        psrplot -p b -x -jT -lpol=0,1 -O -c log=1 -c skip=1 -g 1024x768 -c above:l= -c above:c="Cleaned bandpass (\${type})"     -D \${type}_bandpass.png/png    \$file
    done

    # Create flux and polarisation scrunched archive for SNR images
    pam -Fp -e cleanFp ${cleaned_archive}
    pam -Fp -e rawFp ${raw_archive}
    # Create a frequency, time and polarisation scrunched file for flux calc
    pam -FTp -e cleanFTp ${cleaned_archive}

    # Create matplotlib images and dump the results calculations into a results.json file
    generate_images_results -pid ${project_short} -cleanedfile ${cleaned_archive} -rawfile ${raw_archive} -cleanFp *cleanFp -rawFp *rawFp -cleanFTp *cleanFTp -parfile ${ephemeris} -template ${template} -rcvr ${band} -snr ${snr} -dmfile ${dm_results}
    """
}

process GENERATE_IMAGE_RESULTS_RAW {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/images", mode: 'copy', pattern: "{c,t,r}*png"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/scintillation", mode: 'copy', pattern: "*dynspec*"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: 'copy', pattern: "results.json"

    input:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr)

    output:
    tuple val(pulsar), val(pipe_id), path("*.png")


    """
    # psrplot images
    type=raw
    file=${raw_archive}
    # Do the plots for raw file then cleaned file
    psrplot -p flux -jFTDp -jC                          -g 1024x768 -c above:l= -c above:c="Stokes I Profile (\${type})"     -D \${type}_profile_fts.png/png \$file
    psrplot -p Scyl -jFTD  -jC                          -g 1024x768 -c above:l= -c above:c="Polarisation Profile (\${type})" -D \${type}_profile_ftp.png/png \$file
    psrplot -p freq -jTDp  -jC                          -g 1024x768 -c above:l= -c above:c="Phase vs. Frequency (\${type})"  -D \${type}_phase_freq.png/png  \$file
    psrplot -p time -jFDp  -jC                          -g 1024x768 -c above:l= -c above:c="Phase vs. Time (\${type})"       -D \${type}_phase_time.png/png  \$file
    psrplot -p b -x -jT -lpol=0,1 -O -c log=1 -c skip=1 -g 1024x768 -c above:l= -c above:c="Cleaned bandpass (\${type})"     -D \${type}_bandpass.png/png    \$file

    # Create flux and polarisation scrunched archive for SNR images
    pam -Fp -e rawFp ${raw_archive}

    # Create matplotlib images and dump the results calculations into a results.json file
    generate_images_results -pid ${project_short} -rawfile ${raw_archive} -rawFp *rawFp -parfile ${ephemeris} -rcvr ${band} -snr ${snr}
    """
}


process UPLOAD_RESULTS {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(pulsar), val(pipe_id), path(dat_files), path(png_files), path(dynspec_files), path(results_json)

    """
    #!/usr/bin/env python

    import json
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, decode_id, get_graphql_id
    from psrdb.tables.pipeline_image import PipelineImage
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.toa import Toa

    logger = setup_logging(console=True, level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    pipeline_image_client = PipelineImage(client)
    toa_client            = Toa(client)
    pipeline_run_client   = PipelineRun(client)
    pipeline_run_client.get_dicts = True

    image_data = []
    # grab toa files
    for toa_file in glob("toa*png"):
        if "dmcorrected" in toa_file:
            type = "toa-dm-corrected"
        else:
            type = "toa-single"
        # file_loc, file_type, file_res, cleaned
        image_data.append( (toa_file, type, 'high', True) )

    # file_loc, file_type, file_res, cleaned
    image_data.append( (    "raw_profile_ftp.png",    'profile',     'high', False) )
    image_data.append( ("cleaned_profile_ftp.png",    'profile',     'high', True ) )
    image_data.append( (    "raw_profile_fts.png",    'profile-pol', 'high', False) )
    image_data.append( ("cleaned_profile_fts.png",    'profile-pol', 'high', True ) )
    image_data.append( (    "raw_phase_time.png",     'phase-time',  'high', False) )
    image_data.append( ("cleaned_phase_time.png",     'phase-time',  'high', True ) )
    image_data.append( (    "raw_phase_freq.png",     'phase-freq',  'high', False) )
    image_data.append( ("cleaned_phase_freq.png",     'phase-freq',  'high', True ) )
    image_data.append( (    "raw_bandpass.png",       'bandpass',    'high', False) )
    image_data.append( ("cleaned_bandpass.png",       'bandpass',    'high', True ) )
    image_data.append( (    "raw_SNR_cumulative.png", 'snr-cumul',   'high', False) )
    image_data.append( ("cleaned_SNR_cumulative.png", 'snr-cumul',   'high', True ) )
    image_data.append( (    "raw_SNR_single.png",     'snr-single',  'high', False) )
    image_data.append( ("cleaned_SNR_single.png",     'snr-single',  'high', True ) )

    # Upload images
    for image_path, image_type, resolution, cleaned in image_data:
        image_response = pipeline_image_client.create(
            ${pipe_id},
            image_path,
            image_type,
            resolution,
            cleaned,
        )
        content = json.loads(image_response.content)
        if image_response.status_code not in (200, 201):
            logger.error("Failed to upload image")
            exit(1)

    # Read in results JSON
    with open("results.json", "r") as f:
        results_dict = json.load(f)
    # Update pipeline run as completed (will update pulsarFoldResult)
    pipeline_run_response = pipeline_run_client.update(
        ${pipe_id},
        "Completed",
        results_dict=results_dict,
    )
    """
}


process UPLOAD_RESULTS_RAW {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(pulsar), val(pipe_id), path(png_files)

    """
    #!/usr/bin/env python

    import json
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, decode_id, get_graphql_id
    from psrdb.tables.pipeline_image import PipelineImage
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.toa import Toa

    logger = setup_logging(console=True, level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger)
    pipeline_image_client = PipelineImage(client)
    pipeline_run_client   = PipelineRun(client)
    pipeline_run_client.get_dicts = True

    image_data = []
    # file_loc, file_type, file_res, cleaned
    image_data.append( ("raw_profile_ftp.png",    'profile',     'high', False) )
    image_data.append( ("raw_profile_fts.png",    'profile-pol', 'high', False) )
    image_data.append( ("raw_phase_time.png",     'phase-time',  'high', False) )
    image_data.append( ("raw_phase_freq.png",     'phase-freq',  'high', False) )
    image_data.append( ("raw_bandpass.png",       'bandpass',    'high', False) )
    image_data.append( ("raw_SNR_cumulative.png", 'snr-cumul',   'high', False) )
    image_data.append( ("raw_SNR_single.png",     'snr-single',  'high', False) )

    # Upload images
    for image_path, image_type, resolution, cleaned in image_data:
        image_response = pipeline_image_client.create(
            ${pipe_id},
            image_path,
            image_type,
            resolution,
            cleaned,
        )
        content = json.loads(image_response.content)
        if image_response.status_code not in (200, 201):
            logger.error("Failed to upload image")
            exit(1)

    # Update pipeline run as completed (will update pulsarFoldResult)
    pipeline_run_response = pipeline_run_client.update(
        ${pipe_id},
        "Completed",
        results_dict={
            "percent_rfi_zapped": None,
            "dm": None,
            "dm_err": None,
            "dm_epoch": None,
            "dm_chi2r": None,
            "dm_tres": None,
            "rm": None,
            "rm_err": None,
            "sn": None,
            "flux": None,
        },
    )
    """
}


workflow GENERATE_RESULTS_IMAGES {
    take:
        files_and_meta // channel

    main:
        // Calculate the DM with tempo2 or pdmp
        DM_RM_CALC( files_and_meta.filter { it[8].baseName != "no_template" } )

        // Other images using matplotlib and psrplot and make a results.json
        GENERATE_IMAGE_RESULTS( DM_RM_CALC.out )
        GENERATE_IMAGE_RESULTS_RAW( files_and_meta.filter { it[8].baseName == "no_template" } )

        // Upload images and results
        if ( params.upload ) {
            UPLOAD_RESULTS( GENERATE_IMAGE_RESULTS.out )
            UPLOAD_RESULTS_RAW( GENERATE_IMAGE_RESULTS_RAW.out )
        }
}