/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMeerpipe.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse inputs

// Convert nchan and npols to lists
nchans = params.nchans.split(',').collect { it.toInteger() }
npols  = params.npols.split(',').collect { it.toInteger() }


process MANIFEST_CONFIG_DUMP {
    // Create a json of all the parameters used in this run
    label 'psrdb'

    output:
    path "manifest.json"

    """
    #!/usr/bin/env python

    import json

    manifest = {
        "pipeline_name": "${params.manifest.name}",
        "pipeline_description": "${params.manifest.description}",
        "pipeline_version": "${params.manifest.version}",
        "created_by": "${workflow.userName}",
        "configuration": {
            "utcs": "${params.utcs}",
            "utce": "${params.utce}",
            "obs_pid": "${params.obs_pid}",
            "pulsar": "${params.pulsar}",
            "use_edge_subints": "${params.use_edge_subints}",
            "tos_sn": "${params.tos_sn}",
            "nchans": "${params.nchans}",
            "npols": "${params.npols}",
            "upload": "${params.upload}",
            "psrdb_url": "${params.psrdb_url}",
            "input_dir": "${params.input_dir}",
            "outdir": "${params.outdir}",
            "email": "${params.email}",
            "type": "${params.type}",
            "ephemerides_dir": "${params.ephemerides_dir}",
            "templates_dir": "${params.templates_dir}",
            "ephemeris": "${params.ephemeris}",
            "template": "${params.template}",
        },
    }

    with open("manifest.json", "w") as out_file:
        json.dump(manifest, out_file, indent=4)

    """
}


process OBS_LIST {
    label 'psrdb'
    publishDir "./", mode: 'copy', enabled: params.list_out

    input:
    val utcs
    val utce
    val pulsar
    val obs_pid
    path manifest

    output:
    path "processing_jobs.csv"

    """
    #!/usr/bin/env python

    import os
    import json
    import base64
    import logging
    from datetime import datetime
    from psrdb.tables.observation import Observation
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.template import Template
    from psrdb.tables.ephemeris import Ephemeris
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_rest_api_id, get_graphql_id

    # PSRDB setup
    logger = setup_logging(level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", False, logger=logger)
    obs_client       = Observation(client, "${params.psrdb_token}")
    pipe_run_client  = PipelineRun(client, "${params.psrdb_token}")
    template_client  = Template(   client, "${params.psrdb_token}")
    ephemeris_client = Ephemeris(  client, "${params.psrdb_token}")
    obs_client.get_dicts = True
    obs_client.set_use_pagination(True)

    # Query based on provided parameters
    obs_data = obs_client.list(
        pulsar_name="${pulsar}",
        project_short="${obs_pid}",
        utcs="${utcs}",
        utce="${utce}",
        obs_type="fold",
    )

    # Output file
    with open("processing_jobs.csv", "w") as out_file:
        for ob in obs_data:
            print(ob)
            # Extract data from obs_data
            pulsar   = ob['pulsar']['name']
            obs_id   = int(base64.b64decode(ob['id']).decode("utf-8").split(":")[1])
            utc_obs  = datetime.strptime(ob['utcStart'], '%Y-%m-%dT%H:%M:%S+00:00')
            utc_obs  = "%s-%s" % (utc_obs.date(), utc_obs.time())
            pid_obs  = ob['project']['short']
            pid_code = ob['project']['code']
            beam     = ob['beam']
            band     = ob['band']
            duration = ob['duration']
            cal_loc  = ob['calibration']['location']

            # Grab ephermis and templates
            if "${params.ephemeris}" == "null":
                ephemeris = f"${params.ephemerides_dir}/{pid_obs}/{pulsar}.par"
                if not os.path.exists(ephemeris):
                    # Default to using PTA ephemeris if one does not exist
                    ephemeris = f"${params.ephemerides_dir}/PTA/{pulsar}.par"
            else:
                ephemeris = "${params.ephemeris}"
            if "${params.template}" == "null":
                template = f"${params.templates_dir}/{pid_obs}/{band}/{pulsar}.std"
                if not os.path.exists(template):
                    # Try LBAND template
                    template = f"${params.templates_dir}/{pid_obs}/LBAND/{pulsar}.std"
                if not os.path.exists(template):
                    # Default to using PTA template if one does not exist
                    template = f"${params.templates_dir}/PTA/{band}/{pulsar}.std"
                if not os.path.exists(template):
                    # Final attempt is PTA LBAND template
                    template = f"${params.templates_dir}/PTA/LBAND/{pulsar}.std"
            else:
                template = "${params.template}"

            logger.info(f"Setting up ID: {obs_id} pulsar: {pulsar} band: {band} template: {template} ephemeris: {ephemeris}")

            # Set job as running
            if "${params.upload}" == "true":
                # Get or create template
                template_project = template.split("/")[-3]
                template_band = template.split("/")[-2]
                template_response = template_client.create(
                    pulsar,
                    template_band,
                    template,
                    project_short=template_project,
                )
                logger.debug(template_response)
                template_id = get_rest_api_id(template_response, logging.getLogger(__name__))
                # Get or create ephemeris
                ephemeris_project = ephemeris.split("/")[-2]
                ephemeris_response = ephemeris_client.create(
                    pulsar,
                    ephemeris,
                    project_short=ephemeris_project,
                    comment="",
                )
                logger.debug(ephemeris_response)
                ephemeris_id = get_graphql_id(ephemeris_response, "ephemeris", logging.getLogger(__name__))

                with open("${manifest}", 'r') as file:
                    # Load the JSON data
                    pipeline_config = json.load(file)

                pipe_run_data = pipe_run_client.create(
                    obs_id,
                    ephemeris_id,
                    template_id,
                    "${params.manifest.name}",
                    "${params.manifest.description}",
                    "${params.manifest.version}",
                    "Running",
                    "${params.outdir}",
                    pipeline_config,
                )
                pipe_id = get_graphql_id(pipe_run_data, "pipeline_run", logging.getLogger(__name__))
            else:
                # No uploading so don't make a processing item
                pipe_id = None

            # Write out results
            out_file.write(f"{pulsar},{utc_obs},{pid_obs},{beam},{band},{duration},{cal_loc},{pipe_id},{ephemeris},{template}\\n")
    """
}


process PSRADD_CALIBRATE_CLEAN {
    label 'cpu'
    label 'meerpipe'
    // label 'scratch'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: 'copy', pattern: "*.ar"

    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template)

    output:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path("${pulsar}_${utc}_raw.ar"), path("${pulsar}_${utc}_zap.ar"), env(SNR)

    """
    if ${params.use_edge_subints}; then
        # Grab all archives
        archives=\$(ls ${params.input_dir}/${pulsar}/${utc}/${beam}/*/*.ar)
    else
        if [ -z \$(ls ${params.input_dir}/${pulsar}/${utc}/${beam}/*/*.ar | head -n-1 | tail -n+2) ]; then
            # Grab all archives anyway because there are only two
            archives=\$(ls ${params.input_dir}/${pulsar}/${utc}/${beam}/*/*.ar)
        else
            # Grab all archives except for the first and last one
            archives=\$(ls ${params.input_dir}/${pulsar}/${utc}/${beam}/*/*.ar | head -n-1 | tail -n+2)
        fi
    fi

    echo "Combine the archives"
    psradd -E ${ephemeris} -o ${pulsar}_${utc}_raw.raw \${archives}

    echo "Calibrate the polarisation of the archive"
    if [ "${cal_loc}" == "None" ]; then
        # The archives have already be calibrated so just update the headers
        pac -XP -O ./ -e scalP ${pulsar}_${utc}_raw.raw
    else
        # Use the Stokes paramaters files to calibrate the archive
        pac -Q ${cal_loc} -O ./ -e scal ${pulsar}_${utc}_raw.raw
    fi

    echo "Delay correct"
    dlyfix -e ar ${pulsar}_${utc}_raw.scalP

    echo "Update the RM value if available"
    rm_cat=\$(python -c "from meerpipe.data_load import RM_CAT;print(RM_CAT)")
    if grep -q "${pulsar}" \${rm_cat}; then
        rm=\$(grep ${pulsar} \${rm_cat} | tr -s ' ' | cut -d ' ' -f 2)
        echo "Found RM of \${rm} in the private RM catalogue"
    else
        rm=\$(psrcat -c RM ${pulsar} -X -all | tr -s ' ' | cut -d ' ' -f 1)
        if [ "\${rm}" == "*" ]; then
            echo "No RM found in the ATNF catalogue"
            rm=0
        else
            echo "Found RM of \${rm} in the ATNF catalogue"
        fi
    fi
    pam --RM \${rm} -m ${pulsar}_${utc}_raw.ar

    echo "Clean the archive"
    clean_archive.py -a ${pulsar}_${utc}_raw.ar -T ${template} -o ${pulsar}_${utc}_zap.ar

    # Get the signal to noise ratio of the cleaned archive
    SNR=\$(psrstat -j FTp -c snr=pdmp -c snr ${pulsar}_${utc}_zap.ar | cut -d '=' -f 2)

    echo "Flux calibrate"
    # Create a time and polarisation scruchned profile
    pam -Tp -e tp ${pulsar}_${utc}_zap.ar
    fluxcal -psrname ${pulsar} -obsname ${utc} -obsheader ${params.input_dir}/${pulsar}/${utc}/${beam}/*/obs.header -cleanedfile ${pulsar}_${utc}_zap.ar -rawfile ${pulsar}_${utc}_raw.ar -tpfile *tp -parfile ${ephemeris}
    """
}


process DECIMATE {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/decimated", mode: 'copy', pattern: "${pulsar}_${utc}_zap.*.ar"

    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr)

    output:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path("${pulsar}_${utc}_zap.*.ar")

    """
    for nchan in ${nchans.join(' ')}; do
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

process DM_RM_CALC {
    label 'cpu'
    label 'meerpipe'

    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr)

    output:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path("${pulsar}_${utc}_dm_rm_fit.txt")


    // when:
    // Float.valueOf(snr) > 12.0 // If not enough signal to noise causes tempo2 to core dump

    script:
    if ( Float.valueOf(snr) > 20.0 )
        """
        echo -e "\\nCreate a max channel archive\\n----------------------------------"
        # Calculate nchan to get desired TOA S/N and make sure it is a factor of the archive channels
        arnchan=\$(vap -c nchan ${cleaned_archive} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        nchan=\$(python -c "import math; raw_nchan=math.floor( (${snr}/10.) ** 2); print(next((factor for factor in range(\$arnchan + 1, 2, -1) if \$arnchan % factor == 0 and factor < raw_nchan and factor <= 64)))")
        pam --setnchn \${nchan} -S -p -e dmcalc ${cleaned_archive}
        pam --setnchn \${nchan} -S    -e rmcalc ${cleaned_archive}

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
        # Fit for DM
        tempo2 -nofit -fit DM -set START 40000 -set FINISH 99999 -f ${ephemeris}.dm -outpar ${ephemeris}.dmfit dm.tim

        echo -e "\\nFit for RM\\n----------------------------------"
        input_rm=\$(vap -c rm ${pulsar}_${utc}_zap.rmcalc | tail -n 1| tr -s ' ' | cut -d ' ' -f 2)
        rmfit -D -R \$input_rm -m -100,100,2000 ${pulsar}_${utc}_zap.rmcalc -K /PNG > rmfit_output.txt

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
        pdmp -g ${cleaned_archive}.ps/cps ${cleaned_archive}

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


process GENERATE_TOAS {
    label 'cpu'
    label 'psrchive'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/timing", mode: 'copy', pattern: "*.{residual,tim,par,std}"
    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path(decimated_archives)

    output:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path("*.tim"), path("*.residual")

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
    bash /fred/oz005/users/nswainst/code/meerpipe/tempo2_wrapper.sh \${largest_archive} ${ephemeris}
    """
}
    // DM corrected stuff
    //     echo "Correct for DM\n----------------------------------"
    //     dm=\$(grep DM ${dm_results} | cut -d ' ' -f 2)
    //     pam -D \$dm -e ar.dm_corrected \$ar

    //     echo "Generating TOAs for DM corrected archive\n----------------------------------"
    //     pat -jp \$port  -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM \$ar.dm_corrected  > \${ar}.dm_corrected.tim
    // # And largest DM corrected archive
    // bash /fred/oz005/users/nswainst/code/meerpipe/tempo2_wrapper.sh \${largest_archive}.dm_corrected ${ephemeris}


process GENERATE_IMAGE_RESULTS {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/images", mode: 'copy', pattern: "{c,t,r}*png"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}/scintillation", mode: 'copy', pattern: "*dynspec*"
    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: 'copy', pattern: "results.json"

    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path(dm_results)

    output:
    tuple val(pulsar), val(obs_pid), val(pipe_id), path(ephemeris), path("*.png"), path("*.dat"), path("*dynspec"), path("results.json")


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

    # Created flux and polarisation scrunched archive for SNR images
    pam -Fp -e cleanFp ${cleaned_archive}
    pam -Fp -e rawFp ${raw_archive}

    # Create matplotlib images and dump the results calculations into a results.json file
    generate_images_results -pid ${obs_pid} -cleanedfile ${cleaned_archive} -rawfile ${raw_archive} -cleanFp *cleanFp -rawFp *rawFp -parfile ${ephemeris} -template ${template} -rcvr ${band} -snr ${snr} -dmfile ${dm_results}
    """
}


process UPLOAD_RESULTS {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(pulsar), val(obs_pid), val(pipe_id), path(ephemeris), path(dat_files), path(png_files), path(dynspec_files), path(results_json)

    output:
    tuple val(pulsar), val(obs_pid), val(pipe_id), path(ephemeris)


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
    client = GraphQLClient("${params.psrdb_url}", False, logger)
    pipeline_image_client = PipelineImage(client, "${params.psrdb_token}")
    toa_client            = Toa(client,   "${params.psrdb_token}")
    pipeline_run_client   = PipelineRun(client,   "${params.psrdb_token}")
    pipeline_run_client.set_field_names(True, False)
    pipeline_run_client.get_dicts = True
    pid = '${obs_pid.toLowerCase()}'

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


process UPLOAD_TOAS {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(pulsar), val(utc), val(obs_pid), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr), path(toas), path(residuals)

    output:
    tuple val(pulsar), val(obs_pid), val(pipe_id), path(ephemeris)


    """
    #!/usr/bin/env python

    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, decode_id, get_graphql_id
    from psrdb.tables.pipeline_image import PipelineImage
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.toa import Toa

    logger = setup_logging(console=True, level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", False, logger)
    toa_client            = Toa(client,   "${params.psrdb_token}")
    pipeline_run_client   = PipelineRun(client,   "${params.psrdb_token}")
    pipeline_run_client.set_field_names(True, False)
    pipeline_run_client.get_dicts = True
    pid = '${obs_pid.toLowerCase()}'

    # Upload TOAs
    # Grab ephemeris and template ids
    pipeline_run_data = pipeline_run_client.list(
        id=${pipe_id},
    )
    logger.info(pipeline_run_data)
    ephemeris_id = decode_id(pipeline_run_data[0]["ephemeris"]["id"])
    template_id  = decode_id(pipeline_run_data[0]["template"]["id"])
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
                ephemeris_id,
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
    label 'meerpipe'

    maxForks 1

    input:
    tuple val(pulsar), val(obs_pid), val(pipe_id), path(ephemeris)

    """
    # Loop over each of the TOA filters
    for min_or_max_sub in "--minimum_nsubs" "--maximum_nsubs"; do
        for nchan in ${nchans.join(' ')}; do
            echo -e "\\nDownload the toa file and fit the residuals for \${min_or_max_sub#--} \${nchan}\\n--------------------------\\n"
            #
            psrdb toa download ${pulsar} \$min_or_max_sub --nchan \$nchan
            echo -e "\\nGenerating residuals for \${min_or_max_sub#--} \${nchan}\\n--------------------------\\n"
            bash /fred/oz005/users/nswainst/code/meerpipe/tempo2_wrapper.sh *tim ${ephemeris}
            rm *tim
            if [ -f "toa_${pulsar}_\${min_or_max_sub#--}_nchan\${nchan}.tim.residual" ]; then
                echo -e "\\nUpload the residuals\\n--------------------------\\n"
                psrdb residual create ${pulsar} ${ephemeris} ${obs_pid} toa_${pulsar}_\${min_or_max_sub#--}_nchan\${nchan}.tim.residual
            fi
        done
    done
    """
}

// Info required for completion email and summary
def multiqc_report = []

workflow MEERPIPE {
    MANIFEST_CONFIG_DUMP()

    // Use PSRDB to work out which obs to process
    if ( params.list_in ) {
        // Check contents of list_in
        obs_data = Channel.fromPath( params.list_in ).splitCsv()
    }
    else {
        OBS_LIST(
            params.utcs,
            params.utce,
            params.pulsar,
            params.obs_pid,
            MANIFEST_CONFIG_DUMP.out,
        )
        obs_data = OBS_LIST.out.splitCsv()
    }

    // Combine archives,flux calibrate Clean of RFI with MeerGaurd
    PSRADD_CALIBRATE_CLEAN( obs_data )



    // Calculate the DM with tempo2 or pdmp
    DM_RM_CALC( PSRADD_CALIBRATE_CLEAN.out )

    // Other images using matplotlib and psrplot and make a results.json
    GENERATE_IMAGE_RESULTS( DM_RM_CALC.out )

    // Upload images and results
    if ( params.upload ) {
        UPLOAD_RESULTS( GENERATE_IMAGE_RESULTS.out )
    }



    // Decimate into different time and freq chunnks using pam
    DECIMATE( PSRADD_CALIBRATE_CLEAN.out )

    // Generate TOAs
    GENERATE_TOAS( DECIMATE.out )

    if ( params.upload ) {
        UPLOAD_TOAS( GENERATE_TOAS.out )

        // For each pulsar (not each obs), download all toas and fit residuals
        GENERATE_RESIDUALS( UPLOAD_TOAS.out.groupTuple().map { pulsar, obs_pid, pipe_id, ephemeris -> [ pulsar, obs_pid.first(), pipe_id, ephemeris.first() ] } )
    }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
