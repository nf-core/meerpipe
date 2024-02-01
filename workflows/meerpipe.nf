/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

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
            "ephemeris": "${params.ephemeris}",
            "template": "${params.template}",
        },
    }

    with open("manifest.json", "w") as out_file:
        json.dump(manifest, out_file, indent=4)

    """
}



process GRAB_PREVIOUS_ARCHIVE_SNR {
    label 'process_high'
    label 'meerpipe'

    input:
    tuple val(meta), path(ephemeris), path(template)

    output:
    tuple val(meta), path(ephemeris), path(template), path("empty_raw.ar"), path("${meta.pulsar}_${meta.utc}_zap.ar"), env(SNR), env(FLUX)

    """
    touch empty_raw.ar
    ln -s ${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/${meta.pulsar}_${meta.utc}_zap.ar ${meta.pulsar}_${meta.utc}_zap.ar
    pam -FTp -e FTp ${meta.pulsar}_${meta.utc}_zap.ar
    SNR=\$(psrstat -c snr=pdmp -c snr ${meta.pulsar}_${meta.utc}_zap.FTp | cut -d '=' -f 2)
    FLUX=\$(pdv -f ${meta.pulsar}_${meta.utc}_zap.FTp | tail -n 1 | tr -s ' ' | cut -d ' ' -f 7)
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

    if ( params.obs_csv == "" ) {
        obs_csv = ""
    } else {
        obs_csv = Channel.fromPath(params.obs_csv)
    }
    // Use PSRDB to work out which obs to process
    OBS_LIST(
        params.pulsar,
        params.utcs,
        params.utce,
        params.project,
        obs_csv,
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
        GRAB_PREVIOUS_ARCHIVE_SNR( obs_data.filter { it[2].split('/')[-1] != "no_template.std" } )

        files_and_meta = GRAB_PREVIOUS_ARCHIVE_SNR.out
    } else {
        // Combine archives,flux calibrate Clean of RFI with MeerGaurd
        PSRADD_CALIBRATE_CLEAN( obs_data )

        files_and_meta = PSRADD_CALIBRATE_CLEAN.out
    }

    files_and_meta = files_and_meta
        .map {
            meta, ephemeris, template, raw_archive, cleaned_archive, snr, flux ->
            [
                [
                    id: meta.id,
                    pulsar: meta.pulsar,
                    project_short: meta.project_short,
                    utc: meta.utc,
                    beam: meta.beam,
                    band: meta.band,
                    dur: meta.dur,
                    cal_loc: meta.cal_loc,
                    pipe_id: meta.pipe_id,
                    nchans: meta.nchans,
                    npols: meta.npols,
                    // Two new additions
                    snr: snr,
                    flux: flux,
                ],
                ephemeris,
                template,
                raw_archive,
                cleaned_archive,
            ]
        }

    if ( !params.use_prev_ar ) {
        // Calculate the DM with tempo2 or pdmp
        DM_RM_CALC( files_and_meta )

        // Other images using matplotlib and psrplot and make a results.json
        GENERATE_IMAGE_RESULTS( DM_RM_CALC.out )
        // Upload images and results
        if ( params.upload ) {
            UPLOAD_RESULTS( GENERATE_IMAGE_RESULTS.out )
        }

    }

    // Grab all ephemeris and template pairs for each pulsar
    GRAB_ALL_PAIRS(
        obs_data.map {
            meta, ephemeris, template  ->
            [ meta.pulsar ]
        }
        .unique()
        .collect()
        .flatten()
    )
    pulsar_project_ephem_template = GRAB_ALL_PAIRS.out.splitCsv()

    // Decimate into different time and freq chunks using pam
    DECIMATE(
        files_and_meta
            // Filter out observations without templates
            .filter { it[2].baseName != "no_template" }
            // Only send it the paths and vals it needs
            .map {
                meta, ephemeris, template, raw_archive, cleaned_archive ->
                [ meta, cleaned_archive ]
            }
    )

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
                    id, meta, ephem -> [ meta[0], ephem[0] ]
                }
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
    NfcoreTemplate.dump_parameters(workflow, params)
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
