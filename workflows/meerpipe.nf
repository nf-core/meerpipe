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


process PSRADD_CALIBRATE_CLEAN {
    label 'cpu'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${pulsar}/${utc}/${beam}", mode: params.publish_dir_mode, pattern: "${ template.baseName == "no_template" ? "*raw.ar" : "*zap.ar" }"

    input:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template)

    when:
    beam != "0"

    output:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), path("${pulsar}_${utc}_raw.ar"), path("${pulsar}_${utc}_zap.ar"), env(SNR)

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
    if [[ "${cal_loc}" == "" || "${cal_loc}" == "None" ]]; then
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
        if [[ "\${rm}" == "*" || "\${rm}" == "WARNING*" ]]; then
            echo "No RM found in the ATNF catalogue"
            rm=0
        else
            echo "Found RM of \${rm} in the ATNF catalogue"
        fi
    fi
    pam --RM \${rm} -m ${pulsar}_${utc}_raw.ar

    if [ "${template.baseName}" == "no_template" ]; then
        echo "No template provided so not cleaning archive"
        touch ${pulsar}_${utc}_zap.ar
        SNR=None
    else
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
    fi
    """
}

process GRAB_PREVIOUS_ARCHIVE_SNR {
    label 'cpu'
    label 'meerpipe'

    input:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(cal_loc), val(pipe_id), path(ephemeris), path(template)

    when:
    beam != "0"

    output:
    tuple val(pulsar), val(utc), val(project_short), val(beam), val(band), val(dur), val(pipe_id), path(ephemeris), path(template), env(SNR)

    """
    SNR=\$(psrstat -j FTp -c snr=pdmp -c snr ${params.outdir}/${pulsar}/${utc}/${beam}/${pulsar}_${utc}_zap.ar | cut -d '=' -f 2)
    """
}

// Info required for completion email and summary
def multiqc_report = []

include { OBS_LIST } from '../modules/local/obslist'


include { GENERATE_RESULTS_IMAGES } from '../subworkflows/generate_results_images'
include { DECIMATE_TOA_RESIDUALS } from '../subworkflows/decimate_toa_residuals'

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
    obs_data = OBS_LIST.out.splitCsv()

    if ( params.use_prev_ar ) {
        GRAB_PREVIOUS_ARCHIVE_SNR( obs_data )

        files_and_meta = GRAB_PREVIOUS_ARCHIVE_SNR.out.map {
            pulsar, utc, project_short, beam, band, dur, pipe_id, ephemeris, template, snr ->
            [ pulsar, utc, project_short, beam, band, dur, pipe_id, ephemeris, template, "dummy_file.ar", "${params.outdir}/${pulsar}/${utc}/${beam}/${pulsar}_${utc}_zap.ar", snr ]
        }
    } else {
        // Combine archives,flux calibrate Clean of RFI with MeerGaurd
        PSRADD_CALIBRATE_CLEAN( obs_data )

        files_and_meta = PSRADD_CALIBRATE_CLEAN.out

        // Perform the results and imaging subworkflow which does DM, RM and flux desnsity calculations, creates images and uplaods them
        GENERATE_RESULTS_IMAGES( files_and_meta )
    }

    // Perform the timing subworkflow which does decimation, and creates toas and residuals
    DECIMATE_TOA_RESIDUALS(
        files_and_meta
            // Filter out observations without templates
            .filter { it[8].baseName != "no_template" }
            // Only send it the paths and vals it needs
            .map {
                pulsar, utc, project_short, beam, band, dur, pipe_id, ephemeris, template, raw_archive, cleaned_archive, snr ->
                [ pulsar, utc, beam, dur, pipe_id, cleaned_archive, snr ]
            }
    )
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
