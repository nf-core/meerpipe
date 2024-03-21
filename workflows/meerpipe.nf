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


workflow MEERPIPE {
    if ( params.obs_csv == "" ) {
        obs_csv = Channel.fromPath("none_given")
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
    )

    if ( params.use_prev_ar) {
        // Covert csv into a tupe of the meta map and the files
        files_and_meta = OBS_LIST.out.splitCsv()
        .map {
            pulsar, utc, project_short, beam, band, dur, obs_nchan, obs_nbin, cal_loc, pipe_id, ephemeris, template, n_obs, snr, flux, raw_archive, cleaned_archive ->
            [
                [
                    id: "${pulsar}_${utc}_${beam}",
                    pulsar: pulsar,
                    utc: utc,
                    beam: beam,
                    project_short: project_short,
                    band: band,
                    dur: dur,
                    obs_nchan: obs_nchan,
                    obs_nbin: obs_nbin,
                    cal_loc: cal_loc,
                    pipe_id: pipe_id,
                    nchans: params.nchans.split(',').collect { it.toInteger() },
                    npols:params.npols.split(',').collect  { it.toInteger() },
                    n_obs: n_obs,
                    snr: snr,
                    flux: flux,
                ],
                file(ephemeris),
                file(template),
                file(raw_archive),
                file(cleaned_archive),
            ]
        }
    } else {
        // Covert csv into a tupe of the meta map and the files
        obs_data = OBS_LIST.out.splitCsv()
        .map {
            pulsar, utc, project_short, beam, band, dur, obs_nchan, obs_nbin, cal_loc, pipe_id, ephemeris, template, n_obs ->
            [
                [
                    id: "${pulsar}_${utc}_${beam}",
                    pulsar: pulsar,
                    utc: utc,
                    beam: beam,
                    project_short: project_short,
                    band: band,
                    dur: dur,
                    obs_nchan: obs_nchan,
                    obs_nbin: obs_nbin,
                    cal_loc: cal_loc,
                    pipe_id: pipe_id,
                    nchans: params.nchans.split(',').collect { it.toInteger() },
                    npols:params.npols.split(',').collect  { it.toInteger() },
                    n_obs: n_obs,
                ],
                ephemeris,
                template,
            ]
        }

        // Combine archives,flux calibrate Clean of RFI with MeerGaurd
        PSRADD_CALIBRATE_CLEAN( obs_data )

        files_and_meta = PSRADD_CALIBRATE_CLEAN.out
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
                        obs_nchan: meta.obs_nchan,
                        obs_nbin: meta.obs_nbin,
                        cal_loc: meta.cal_loc,
                        pipe_id: meta.pipe_id,
                        nchans: meta.nchans,
                        npols: meta.npols,
                        n_obs: meta.n_obs,
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
        OBS_LIST.out
            .splitCsv()
            .map { [ it[0] ] }
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
                        obs_nchan: meta.obs_nchan,
                        obs_nbin: meta.obs_nbin,
                        cal_loc: meta.cal_loc,
                        pipe_id: meta.pipe_id,
                        nchans: meta.nchans,
                        npols: meta.npols,
                        n_obs: meta.n_obs,
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
                    [ groupKey( "${meta.pulsar}_${meta.project_short}", meta.n_obs.toInteger() ), meta, ephem ]
                }
                .groupTuple( remainder: true )
                .map {
                    id, meta, ephem -> [ meta[0], ephem[0] ]
                }
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
