process GENERATE_IMAGE_RESULTS {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'
    label 'scratch'

    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/images", mode: 'copy', pattern: "{c,t,r}*png"
    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}/scintillation", mode: 'copy', pattern: "*dynspec*"
    publishDir "${params.outdir}/${meta.pulsar}/${meta.utc}/${meta.beam}", mode: 'copy', pattern: "results.json"

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meerpipe:latest':
        'nickswainston/meerpipe:3.0.6' }"

    input:
    tuple val(meta), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), path(dm_results), path(rm_fit_image)

    output:
    tuple val(meta), path("*.png", includeInputs: true), path("*.dat"), path("*dynspec"), path("results.json")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # psrplot images
    raw_archive=${raw_archive}
    cleaned_archive=${cleaned_archive}
    type_file_array=(${ template.baseName == "no_template" ? '"raw ${raw_archive}"' : '"raw ${raw_archive}" "cleaned ${cleaned_archive}"' })
    for i in "\${type_file_array[@]}"; do
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
    pam -Fp -e rawFp ${raw_archive}
    if [ "${template.baseName}" != "no_template" ]; then
        pam -Fp -e cleanFp ${cleaned_archive}
        # Create a frequency, time and polarisation scrunched file for flux calc
        pam -FTp -e cleanFTp ${cleaned_archive}

        echo "Check if you need to change the template bins"
        obs_nbin=\$(vap -c nbin ${cleaned_archive} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        std_nbin=\$(vap -c nbin ${template} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        if [ "\$obs_nbin" == "\$std_nbin" ]; then
            std_template=${template}
        else
            echo "Making a new template with right number of bins"
            pam -b \$((std_nbin / obs_nbin)) -e new_std ${template}
            std_template=*new_std
        fi
        psrflux -s \${std_template} -e dynspec ${raw_archive}
        psrflux -s \${std_template} -e dynspec ${cleaned_archive}
    fi

    # Create matplotlib images and dump the results calculations into a results.json file
    generate_images_results \\
        --pid ${meta.project_short} \\
        --raw_file ${raw_archive} \\
        --raw_Fp *rawFp \\
        --cleaned_file ${cleaned_archive} \\
        --clean_Fp *cleanFp \\
        --clean_FTp *cleanFTp \\
        --par_file ${ephemeris} \\
        --template ${template} \\
        --rcvr ${meta.band} \\
        --snr ${meta.snr} \\
        --flux ${meta.flux} \\
        --dm_file ${dm_results} \\
        ${ template.baseName == "no_template" ? "--raw_only" : "" }
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.png"
    touch "${prefix}.dat"
    touch "${prefix}.dynspec"
    touch results.json
    """
}
