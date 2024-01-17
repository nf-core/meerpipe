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
            "ephemeris": "${params.ephemeris}",
            "template": "${params.template}",
        },
    }

    with open("manifest.json", "w") as out_file:
        json.dump(manifest, out_file, indent=4)

    """
}



process GRAB_PREVIOUS_ARCHIVE_SNR {
    label 'cpu'
    label 'meerpipe'

    input:
    tuple val(meta), path(ephemeris), path(template)

    output:
    tuple val(meta), path(ephemeris), path(template), env(SNR)

    """
    if [ "${template.baseName}" == "no_template" ]; then
        SNR=None
    else
        SNR=\$(psrstat -j FTp -c snr=pdmp -c snr ${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/${meta.pulsar}_${meta.utc}_zap.ar | cut -d '=' -f 2)
    fi
    """
}




process UPLOAD_RESULTS_RAW {
    label 'psrdb'

    maxForks 1

    input:
    tuple val(meta), path(png_files)

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
            ${meta.pipe_id},
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
        ${meta.pipe_id},
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


process GENERATE_IMAGE_RESULTS_RAW {
    label 'cpu'
    label 'meerpipe'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/images", mode: 'copy', pattern: "{c,t,r}*png"
    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/scintillation", mode: 'copy', pattern: "*dynspec*"
    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}", mode: 'copy', pattern: "results.json"

    input:
    tuple val(meta), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), val(snr)

    output:
    tuple val(meta), path("*.png")


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
    generate_images_results -pid ${meta.project_short} -rawfile ${raw_archive} -rawFp *rawFp -parfile ${ephemeris} -rcvr ${meta.band} -snr ${snr}
    """
}

// Info required for completion email and summary
def multiqc_report = []

include { OBS_LIST               } from '../modules/local/obs_list'
include { PSRADD_CALIBRATE_CLEAN } from '../modules/local/psradd_calibrate_clean'
include { DM_RM_CALC             } from '../modules/local/dm_rm_calc'
include { GENERATE_IMAGE_RESULTS } from '../modules/local/generate_image_results'
include { UPLOAD_RESULTS         } from '../modules/local/upload_results'
include { GRAB_ALL_PAIRS         } from '../modules/local/grab_all_pairs'
include { DECIMATE               } from '../modules/local/decimate'
include { GENERATE_TOAS          } from '../modules/local/generate_toas'
include { UPLOAD_TOAS            } from '../modules/local/upload_toas'
include { GENERATE_RESIDUALS     } from '../modules/local/generate_residuals'


workflow MEERPIPE {
    MANIFEST_CONFIG_DUMP()

    // Use PSRDB to work out which obs to process
    OBS_LIST(
        params.pulsar,
        params.utcs,
        params.utce,
        params.project,
        params.obs_csv,
        params.upload,
        params.psrdb_url,
        params.psrdb_token,
        params.ephemeris,
        params.template,
        params.outdir,
        MANIFEST_CONFIG_DUMP.out,
    )

    // Covert csv into a tupe of the meta map and the files
    obs_data = OBS_LIST.out.splitCsv()
    .map {
        pulsar, utc, project_short, beam, band, dur, cal_loc, pipe_id, ephemeris, template ->
        [
            [
                id: "${pulsar}_${utc}_${beam}",
                pulsar: pulsar,
                utc: utc,
                beam: beam,
                project_short: project_short,
                band: band,
                dur: dur,
                cal_loc: cal_loc,
                pipe_id: pipe_id,
                nchans: params.nchans.split(',').collect { it.toInteger() },
                npols:params.npols.split(',').collect  { it.toInteger() },
            ],
            ephemeris,
            template,
        ]
    }

    if ( params.use_prev_ar ) {
        GRAB_PREVIOUS_ARCHIVE_SNR( obs_data )

        files_and_meta = GRAB_PREVIOUS_ARCHIVE_SNR.out.map {
            pulsar, utc, project_short, beam, band, dur, pipe_id, ephemeris, template, snr ->
            [ pulsar, utc, project_short, beam, band, dur, pipe_id, ephemeris, template, "dummy_file.ar", "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/${meta.pulsar}_${meta.utc}_zap.ar", snr ]
        }
    } else {
        // Combine archives,flux calibrate Clean of RFI with MeerGaurd
        PSRADD_CALIBRATE_CLEAN( obs_data )

        files_and_meta = PSRADD_CALIBRATE_CLEAN.out
    }

    // Calculate the DM with tempo2 or pdmp
    DM_RM_CALC( files_and_meta.filter { it[2].baseName != "no_template" } )

    // Other images using matplotlib and psrplot and make a results.json
    GENERATE_IMAGE_RESULTS( DM_RM_CALC.out )
    GENERATE_IMAGE_RESULTS_RAW( files_and_meta.filter { it[2].baseName == "no_template" } )

    // Upload images and results
    if ( params.upload ) {
        UPLOAD_RESULTS( GENERATE_IMAGE_RESULTS.out )
        UPLOAD_RESULTS_RAW( GENERATE_IMAGE_RESULTS_RAW.out )
    }

    // Perform the timing subworkflow which does decimation, and creates toas and residuals
    dtr_files_and_meta = files_and_meta
        // Filter out observations without templates
        .filter { it[2].baseName != "no_template" }
        // Only send it the paths and vals it needs
        .map {
            meta, ephemeris, template, raw_archive, cleaned_archive, snr ->
            [ meta, cleaned_archive, snr ]
        }
    // Grab all ephemeris and template pairs for each pulsar
    GRAB_ALL_PAIRS(
        dtr_files_and_meta.map {
            meta, cleaned_archive, snr  ->
            [ meta.pulsar ]
        }
        .unique()
        .collect()
        .flatten()
    )
    pulsar_project_ephem_template = GRAB_ALL_PAIRS.out.splitCsv()

    // Decimate into different time and freq chunks using pam
    DECIMATE( dtr_files_and_meta )

    // Generate TOAs
    GENERATE_TOAS(
        pulsar_project_ephem_template
            .combine(
                DECIMATE
                .out
                .map {
                    meta, decimated_archives ->
                    [ meta.pulsar, meta, decimated_archives ]
                },
                by: 0
            )
            .map {
                pulsar, project_short, ephemeris, template, meta, decimated_archives ->
                [
                    [
                        id: "${pulsar}_${meta.utc}_${meta.beam}_${project_short}",
                        pulsar: pulsar,
                        project_short: project_short,
                        utc: meta.utc,
                        beam: meta.beam,
                        band: meta.band,
                        dur: meta.dur,
                        cal_loc: meta.cal_loc,
                        pipe_id: meta.pipe_id,
                        nchans: meta.nchans,
                        npols: meta.npols,
                    ],
                    ephemeris,
                    template,
                    decimated_archives,
                ]
            }
    )

    if ( params.upload ) {
        UPLOAD_TOAS( GENERATE_TOAS.out )

        // For each pulsar (not each obs), download all toas and fit residuals
        GENERATE_RESIDUALS(
            UPLOAD_TOAS
                .out
                .map {
                    meta, ephem ->
                    [ "${meta.pulsar}_${meta.project_short}", meta, ephem ]
                }
                .groupTuple()
                .map {
                    id, meta, ephem -> [ meta[0], ephem ]
                }.view()
        )
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
