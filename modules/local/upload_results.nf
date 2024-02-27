// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process UPLOAD_RESULTS {
    tag "$meta.id"
    label 'psrdb'
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(png_files), path(dat_files), path(dynspec_files), path(results_json)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #!/usr/bin/env python

    import json
    import time
    import logging
    from glob import glob
    from psrdb.graphql_client import GraphQLClient
    from psrdb.utils.other import setup_logging, decode_id, get_graphql_id
    from psrdb.tables.pipeline_image import PipelineImage
    from psrdb.tables.pipeline_run import PipelineRun
    from psrdb.tables.toa import Toa

    logger = setup_logging(console=True, level=logging.DEBUG)
    if ${task.attempt} > 1:
        wait_time = (${task.attempt} - 1) * 30
        logger.info(f"Waiting for \${wait_time} s before retrying")
        time.sleep(wait_time)

    raw_only = ${ dynspec_files.baseName == "empty" ? "True" : "False" }

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
    image_data.append( ("raw_profile_ftp.png",    'profile',     'high', False) )
    image_data.append( ("raw_profile_fts.png",    'profile-pol', 'high', False) )
    image_data.append( ("raw_phase_time.png",     'phase-time',  'high', False) )
    image_data.append( ("raw_phase_freq.png",     'phase-freq',  'high', False) )
    image_data.append( ("raw_bandpass.png",       'bandpass',    'high', False) )
    image_data.append( ("raw_SNR_cumulative.png", 'snr-cumul',   'high', False) )
    image_data.append( ("raw_SNR_single.png",     'snr-single',  'high', False) )
    image_data.append( ("${meta.pulsar}_${meta.utc}_raw.ar.dynspec.png", 'dynamic-spectrum', 'high', False) )
    if not raw_only:
        image_data.append( ("cleaned_profile_ftp.png",    'profile',     'high', True ) )
        image_data.append( ("cleaned_profile_fts.png",    'profile-pol', 'high', True ) )
        image_data.append( ("cleaned_phase_time.png",     'phase-time',  'high', True ) )
        image_data.append( ("cleaned_phase_freq.png",     'phase-freq',  'high', True ) )
        image_data.append( ("cleaned_bandpass.png",       'bandpass',    'high', True ) )
        image_data.append( ("cleaned_SNR_cumulative.png", 'snr-cumul',   'high', True ) )
        image_data.append( ("${meta.pulsar}_${meta.utc}_zap.ar.dynspec.png", 'dynamic-spectrum', 'high', True ) )

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

    # Read in results JSON
    with open("results.json", "r") as f:
        results_dict = json.load(f)
    # Update pipeline run as completed (will update pulsarFoldResult)
    pipeline_run_response = pipeline_run_client.update(
        ${meta.pipe_id},
        "Completed",
        results_dict=results_dict,
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch done
    """
}
