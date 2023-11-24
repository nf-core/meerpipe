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
npols  = params.npols.split(',').collect  { it.toInteger() }


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
            "project": "${params.project}",
            "obs_csv": "${params.obs_csv}",
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

    input:
    val utcs
    val utce
    val pulsars
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
    import pandas as pd
    from datetime import datetime
    from psrdb.tables.observation import Observation
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.template import Template
    from psrdb.tables.ephemeris import Ephemeris
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, get_rest_api_id, get_graphql_id, decode_id

    # PSRDB setup
    logger = setup_logging(level=logging.DEBUG)
    client = GraphQLClient("${params.psrdb_url}", "${params.psrdb_token}", logger=logger)
    obs_client       = Observation(client)
    pipe_run_client  = PipelineRun(client)
    template_client  = Template(client)
    ephemeris_client = Ephemeris(client)
    obs_client.get_dicts = True
    obs_client.set_use_pagination(True)

    if "${params.obs_csv}" == "null":
        # Query based on provided parameters
        obs_data = obs_client.list(
            pulsar_name=["${pulsars.split(',').join('","')}"],
            project_short="${obs_pid}",
            utcs="${utcs}",
            utce="${utce}",
            obs_type="fold",
        )
        obs_df = pd.DataFrame(columns=["Obs ID","Pulsar Jname","UTC Start","Project Short Name","Beam #","Observing Band","Duration (s)","Calibration Location"])
        for obs in obs_data:
            obs_df = pd.concat(
                [
                    obs_df,
                    pd.Series({
                        "Obs ID": decode_id(obs['id']),
                        "Pulsar Jname": obs['pulsar']['name'],
                        "UTC Start": datetime.strptime(
                            obs['utcStart'],
                            '%Y-%m-%dT%H:%M:%S+00:00',
                        ).strftime('%Y-%m-%d-%H:%M:%S'),
                        "Project Short Name": obs['project']['short'],
                        "Beam #": obs['beam'],
                        "Observing Band": obs['band'],
                        "Duration (s)": obs['duration'],
                        "Calibration Location": obs['calibration']['location'],
                    }).to_frame().T,
                ],
                ignore_index=True
            )
    else:
        # Read in obs from csv
        obs_df = pd.read_csv("${params.obs_csv}")


    # Grab all the ephems and templates for each pulsar first
    pulsar_ephem_template = {}
    for _, obs in obs_df.iterrows():
        # Extract data from obs_df
        pulsar   = obs['Pulsar Jname']
        project  = obs['Project Short Name']
        band     = obs['Observing Band']
        if pulsar in pulsar_ephem_template.keys():
            if band in pulsar_ephem_template[pulsar].keys():
                # Already grabbed so continue
                continue
        else:
            pulsar_ephem_template[pulsar] = {}

        # Grab ephermis and templates
        if "${params.ephemeris}" == "null":
            ephemeris = f"${params.ephemerides_dir}/{project}/{pulsar}.par"
            if not os.path.exists(ephemeris):
                # Default to using PTA ephemeris if one does not exist
                ephemeris = f"${params.ephemerides_dir}/PTA/{pulsar}.par"
            if not os.path.exists(ephemeris):
                # Default to using TPA ephemeris if one does not exist
                ephemeris = f"${params.ephemerides_dir}/TPA/{pulsar}.par"
            if not os.path.exists(ephemeris):
                # Default to using RelBin ephemeris if one does not exist
                ephemeris = f"${params.ephemerides_dir}/RelBin/{pulsar}.par"
            if not os.path.exists(ephemeris):
                # Default to using PTUSE ephemeris if one does not exist
                ephemeris = f"${params.ephemerides_dir}/PTUSE/{pulsar}.par"
        else:
            ephemeris = "${params.ephemeris}"
        if "${params.template}" == "null":
            template = f"${params.templates_dir}/{project}/{band}/{pulsar}.std"
            if not os.path.exists(template):
                # Try LBAND template
                template = f"${params.templates_dir}/{project}/LBAND/{pulsar}.std"
            if not os.path.exists(template):
                # Default to using PTA template if one does not exist
                template = f"${params.templates_dir}/PTA/{band}/{pulsar}.std"
            if not os.path.exists(template):
                # Final attempt is PTA LBAND template
                template = f"${params.templates_dir}/PTA/LBAND/{pulsar}.std"
            if not os.path.exists(template):
                # Default to using TPA template if one does not exist
                template = f"${params.templates_dir}/TPA/{band}/{pulsar}.std"
            if not os.path.exists(template):
                # Final attempt is TPA LBAND template
                template = f"${params.templates_dir}/TPA/LBAND/{pulsar}.std"
            if not os.path.exists(template):
                # Default to using RelBin template if one does not exist
                template = f"${params.templates_dir}/RelBin/{band}/{pulsar}.std"
            if not os.path.exists(template):
                # Final attempt is RelBin LBAND template
                template = f"${params.templates_dir}/RelBin/LBAND/{pulsar}.std"
        else:
            template = "${params.template}"
        ephem_template = {
            "ephemeris": ephemeris,
            "template": template,
        }
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
            ephem_template["template_id"] = get_rest_api_id(template_response, logging.getLogger(__name__))
            # Get or create ephemeris
            ephemeris_project = ephemeris.split("/")[-2]
            ephemeris_response = ephemeris_client.create(
                pulsar,
                ephemeris,
                project_short=ephemeris_project,
                comment="",
            )
            logger.debug(ephemeris_response)
            ephem_template["ephemeris_id"] = get_graphql_id(ephemeris_response, "ephemeris", logging.getLogger(__name__))
        pulsar_ephem_template[pulsar][band] = ephem_template



    # Add ephemeris and template to df and create the pipeline run object
    obs_df['pipe_id'] = ''
    obs_df['ephemeris'] = ''
    obs_df['template'] = ''
    for index, obs in obs_df.iterrows():
        ephemeris = pulsar_ephem_template[pulsar][band]["ephemeris"]
        template  = pulsar_ephem_template[pulsar][band]["template"]

        logger.info(f"Setting up ID: {obs['Obs ID']} pulsar: {obs['Pulsar Jname']} band: {obs['Observing Band']} template: {template} ephemeris: {ephemeris}")

        # Set job as running
        if "${params.upload}" == "true":
            ephemeris_id = pulsar_ephem_template[pulsar][band]["ephemeris_id"]
            template_id  = pulsar_ephem_template[pulsar][band]["template_id"]
            with open("${manifest}", 'r') as file:
                # Load the JSON data
                pipeline_config = json.load(file)

            pipe_run_data = pipe_run_client.create(
                obs['Obs ID'],
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
        obs_df.at[index, 'pipe_id']   = pipe_id
        obs_df.at[index, 'ephemeris'] = ephemeris
        obs_df.at[index, 'template']  = template

    # Write out results
    obs_df.drop('Obs ID', axis=1, inplace=True)
    obs_df.to_csv("processing_jobs.csv", header=False, index=False)
    """
}


process PSRADD_CALIBRATE_CLEAN {
    label 'cpu'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: 'copy', pattern: "*zap.ar"

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
    if [ "${cal_loc}" == "" ]; then
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

    echo "Check if you need to change the template bins"
    obs_nbin=\$(vap -c nbin ${pulsar}_${utc}_raw.raw | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
    std_nbin=\$(vap -c nbin ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
    if [ "\$obs_nbin" == "\$std_nbin"]; then
        std_template=${template}
    else
        echo "Making a new template with right number of bins"
        pam -b \$((std_nbin / obs_nbin)) -e new_std ${template}
        std_template=*new_std
    fi
    echo "Clean the archive"
    clean_archive.py -a ${pulsar}_${utc}_raw.ar -T \${std_template} -o ${pulsar}_${utc}_zap.ar

    # Get the signal to noise ratio of the cleaned archive
    SNR=\$(psrstat -j FTp -c snr=pdmp -c snr ${pulsar}_${utc}_zap.ar | cut -d '=' -f 2)

    echo "Flux calibrate"
    # Create a time and polarisation scruchned profile
    pam -Tp -e tp ${pulsar}_${utc}_zap.ar
    fluxcal -psrname ${pulsar} -obsname ${utc} -obsheader ${params.input_dir}/${pulsar}/${utc}/${beam}/*/obs.header -cleanedfile ${pulsar}_${utc}_zap.ar -rawfile ${pulsar}_${utc}_raw.ar -tpfile *tp -parfile ${ephemeris}
    """
}

// Info required for completion email and summary
def multiqc_report = []

include { GENERATE_RESULTS_IMAGES } from '../subworkflows/generate_results_images'
include { DECIMATE_TOA_RESIDUALS } from '../subworkflows/decimate_toa_residuals'

workflow MEERPIPE {
    MANIFEST_CONFIG_DUMP()

    // Use PSRDB to work out which obs to process
    OBS_LIST(
        params.utcs,
        params.utce,
        params.pulsar,
        params.project,
        MANIFEST_CONFIG_DUMP.out,
    )
    obs_data = OBS_LIST.out.splitCsv()

    // Combine archives,flux calibrate Clean of RFI with MeerGaurd
    PSRADD_CALIBRATE_CLEAN( obs_data )

    GENERATE_RESULTS_IMAGES( PSRADD_CALIBRATE_CLEAN.out )

    DECIMATE_TOA_RESIDUALS( PSRADD_CALIBRATE_CLEAN.out )
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
