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

process DM_RM_CALC {
    tag "$meta.id"
    label 'process_high'
    label 'meerpipe'

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
    tuple val(meta), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive)

    output:
    tuple val(meta), path(ephemeris), path(template), path(raw_archive), path(cleaned_archive), path("${meta.pulsar}_${meta.utc}_dm_rm_fit.json")

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
     if ( task.attempt > 2 || template.baseName == "no_template" )
        """
        cat << EOF > ${meta.pulsar}_${meta.utc}_dm_rm_fit.json
        {
            "DM": "None",
            "ERR": "None",
            "EPOCH": "None",
            "CHI2R": "None",
            "TRES": "None",
            "RM": "None",
            "RM_ERR": "None"
        }
        EOF
        """
    else if ( Float.valueOf(meta.snr) > 20.0 )
        """
        echo -e "\\nCreate a max channel archive\\n----------------------------------"
        # Calculate nchan to get desired TOA S/N and make sure it is a factor of the archive channels
        arnchan=\$(vap -c nchan ${cleaned_archive} | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2)
        nchan=\$(python -c "import math; raw_nchan=math.floor( (${meta.snr}/10.) ** 2); print(next((factor for factor in range(\$arnchan + 1, 2, -1) if \$arnchan % factor == 0 and factor < raw_nchan and factor <= 64)))")
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
        pat -jp \$port -f "tempo2 IPTA" -C "chan rcvr snr length subint" -s ${template} -A FDM ${meta.pulsar}_${meta.utc}_zap.dmcalc > dm.tim

        echo -e "\\nCalc DM with tempo2\\n----------------------------------"
        # Remove dm derivatives
        sed '/^DM[1-9]/d' ${ephemeris} > ${ephemeris}.dm
        echo "MODE 1" >>  ${ephemeris}.dm
        # Remove zero S/N TOAs
        sed -i '/-snr 0 /d' dm.tim
        # Fit for DM
        tempo2 -nofit -fit DM -set START 40000 -set FINISH 99999 -f ${ephemeris}.dm -outpar ${ephemeris}.dmfit dm.tim

        input_rm=\$(vap -c rm ${meta.pulsar}_${meta.utc}_zap.rmcalc | tail -n 1| tr -s ' ' | cut -d ' ' -f 2)
        lower_rm=\$(echo "\$input_rm - 34" | bc -l)
        higher_rm=\$(echo "\$input_rm + 34" | bc -l)
        echo -e "\\nFit for RM from \$lower_rm - \$higher_rm \\n----------------------------------"
        rmfit -D -m \$lower_rm,\$higher_rm,200 ${meta.pulsar}_${meta.utc}_zap.rmcalc -K /PNG > rmfit_output.txt

        echo -e "\\nGrab the outputs and write it to a file\\n----------------------------------"
        DM=\$(grep    "^DM "      ${ephemeris}.dmfit | awk '{print \$2}')
        ERR=\$(grep   "^DM "      ${ephemeris}.dmfit | awk '{print \$4}')
        EPOCH=\$(grep "^DMEPOCH " ${ephemeris}.dmfit | awk '{print \$2}')
        CHI2R=\$(grep "^CHI2R "   ${ephemeris}.dmfit | awk '{print \$2}')
        TRES=\$(grep  "^TRES "    ${ephemeris}.dmfit | awk '{print \$2}')
        rm_results=\$(grep "Best RM is" rmfit_output.txt | cut -d ':' -f 2)
        RM=\$(echo     \$rm_results | cut -d '/' -f 1 | cut -d ' ' -f 1)
        if grep -q "+/-" "rmfit_output.txt"; then
            RM_ERR=\$(echo \$rm_results | cut -d '/' -f 2 | cut -d ' ' -f 2)
        else
            echo "rmfit do not return an uncertainty so recording it as None"
            RM_ERR=None
        fi

        cat << EOF > ${meta.pulsar}_${meta.utc}_dm_rm_fit.json
        {
            "DM": "\${DM}",
            "ERR": "\${ERR}",
            "EPOCH": "\${EPOCH}",
            "CHI2R": "\${CHI2R}",
            "TRES": "\${TRES}",
            "RM": "\${RM}",
            "RM_ERR": "\${RM_ERR}"
        }
        EOF
        """
    else
        """
        pam -T -S -p -e dmcalc ${cleaned_archive}
        pdmp -g ${cleaned_archive}.ps/cps -mc 16 *dmcalc

        # Grab the outputs and write it to a file
        DM=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 5)
        ERR=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 6)
        EPOCH=\$(cat pdmp.per | tr -s ' ' | cut -d ' ' -f 2)

        cat << EOF > ${meta.pulsar}_${meta.utc}_dm_rm_fit.json
        {
            "DM": "\${DM}",
            "ERR": "\${ERR}",
            "EPOCH": "\${EPOCH}",
            "CHI2R": "None",
            "TRES": "None",
            "RM": "None",
            "RM_ERR": "None"
        }
        EOF
        """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    cat << EOF > ${meta.pulsar}_${meta.utc}_dm_rm_fit.json
    {
        "DM": "None",
        "ERR": "None",
        "EPOCH": "None",
        "CHI2R": "None",
        "TRES": "None",
        "RM": "None",
        "RM_ERR": "None"
    }
    EOF
    """
}
