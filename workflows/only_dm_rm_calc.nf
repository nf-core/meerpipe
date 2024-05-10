/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OBS_LIST             } from '../modules/local/obs_list'
include { DM_RM_CALC           } from '../modules/local/dm_rm_calc'
include { UPLOAD_DM_RM_RESULTS } from '../modules/local/upload_dm_rm_results'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ONLY_DM_RM_CALC {
    if ( params.obs_csv == "" ) {
        obs_csv = Channel.fromPath("none_given")
    } else {
        obs_csv = Channel.fromPath(params.obs_csv)
    }

    if ( !params.use_prev_ar) {
        println("Must use --use_prev_ar for this pipeline")
        exit(1)
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
    )

    // Covert csv into a tupe of the meta map and the files
    files_and_meta = OBS_LIST.out.splitCsv()
    .map {
        pulsar, utc, project_short, beam, band, dur, mode_dur, obs_nchan, obs_nbin, cal_loc, pipe_id, ephemeris, template, n_obs, snr, flux, percent_rfi_zapped, raw_archive, cleaned_archive ->
        [
            [
                id: "${pulsar}_${utc}_${beam}",
                pulsar: pulsar,
                utc: utc,
                beam: beam,
                project_short: project_short,
                band: band,
                dur: dur,
                mode_dur: mode_dur,
                obs_nchan: obs_nchan,
                obs_nbin: obs_nbin,
                pipe_id: pipe_id,
                nchans: params.nchans.split(',').collect { it.toInteger() },
                npols:params.npols.split(',').collect  { it.toInteger() },
                n_obs: n_obs,
                snr: snr,
                flux: flux,
                percent_rfi_zapped: percent_rfi_zapped,
            ],
            file(ephemeris),
            file(template),
            file(raw_archive),
            file(cleaned_archive),
        ]
    }

    // Calculate the DM with tempo2 or pdmp
    DM_RM_CALC( files_and_meta )

    // Upload images and results
    if ( params.upload ) {
        UPLOAD_DM_RM_RESULTS(
            DM_RM_CALC.out.map {
                meta, ephemeris, template, raw_archive, cleaned_archive, results_json, rm_image ->
                [ meta, results_json, rm_image ]
            }
        )
    }
}

workflow {
    ONLY_DM_RM_CALC()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
