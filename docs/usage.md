# nf-core/meerpipe: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/meerpipe/usage](https://nf-co.re/meerpipe/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

A processing pipeline to convert raw pulsar archives from the MeerKAT telescopes into timing data and other useful products.

## Archive input

The pipeline works out what observations are available based on the input filters and then find the files in the directory structure.
The directory structure is defined as:

```console
${params.input_dir}/${meta.pulsar}/${meta.utc}/${meta.beam}/*/*.ar
```

This structure is based on what is currently available on the Swinburne's OzSTAR supercomputer.

```console
/fred/oz005/timing/<pulsar>/<utc>/<beam>/<freq>
```

where:

 - `<pulsar>`: A pulsar J name (e.g. J0437-4715)
 - `<utc>`: The start time of the observation in UTC and the format "YYYY-MM-DD-HH:MM:SS.SS" (e.g. 2023-12-11-03:23:30)
 - `<beam>`: The beam ID (the PTUSE server, e.g. 4)
 - `<freq>`: The centre frequency in MHz (e.g. 1284)

 These observations you wish to process can be defined using either the arguments on the command line ([Filtering observations](/usage#filtering-observations)) or with a observation file for more advanced combinations of observations ([Observation file](usage#observation-file)).


### Filtering observations

The pipeline can filter observations based on the [Observation selection](/parameters#observation-selection) parameters.
Some common ways to use the filters include processing all observations for a single pulsar:

```console
--pulsar <pulsar>
```

or to process a single observation, you can use the same date value for `--utcs` and `--utce`:

```console
--pulsar <pulsar> --utcs <utc> --utce <utc>
```

### Observation file

The observation file is a CSV file with the following columns:

| Column    | Description |
| --------- | ----------- |
| Obs ID | The observation's database ID |
| Pulsar Jname | The pulsar's Jname |
| UTC Start | The start time of the observation at UTC in the format "YYYY-MM-DD-HH:MM:SS.SS"  |
| Project Short Name | The abbreviated name of the observing project (e.g. PTA, TPA, Relbin or GC) |
| Beam # | An integer from 1 to 4 of the PTUSE machine ID that preprocessed the observations|
| Observing Band | The frequency band of the observation (e.g UHF, LBAND or SBAND_<0-4>) |
| Duration (s) | Duration of the observation in seconds |
| Calibration Location | The path to the calibration file if it exists |

For example, one of the observation files used for testing looks like this:

```console
Obs ID,Pulsar Jname,UTC Start,Project Short Name,Beam #,Observing Band,Duration (s),Calibration Location
35409,J1534-5334,2023-05-05-01:07:10,TPA,1,LBAND,116.53766938317756,None
9017,J1418-3921,2020-08-08-11:55:30,TPA,4,LBAND,233.7823401869157,None
2954,J1013-5934,2020-01-04-20:29:13,TPA,1,LBAND,457.50813667289714,/fred/oz005/users/aparthas/reprocessing_MK/poln_calibration/2020-01-04-18:56:54.jones
10173,J1410-7404,2020-09-05-10:00:47,TPA,2,LBAND,89.60730437383191,None
173,J0437-4715,2019-03-26-16:26:02,PTA,1,LBAND,15.999999999999988,/fred/oz005/users/aparthas/reprocessing_MK/poln_calibration/2019-03-26-16:10:39.jones
12929,J1919+0021,2020-11-30-14:55:13,TPA,4,LBAND,1318.9994574953269,None
3608,J0514-4408,2020-02-21-17:50:16,TPA,1,LBAND,1799.9996901682246,/fred/oz005/users/aparthas/reprocessing_MK/poln_calibration/2020-02-21-17:32:22.jones
16650,J0900-3144,2021-04-20-14:47:54,RelBin,3,UHF,2047.541722352941,None
1314,J1811-2405,2019-09-24-13:51:57,RelBin,1,LBAND,7205.578322542057,/fred/oz005/users/aparthas/reprocessing_MK/poln_calibration/2019-09-24-11:31:01.jones
1558,J0955-6150,2019-10-07-02:26:56,RelBin,1,LBAND,14402.595572037384,/fred/oz005/users/aparthas/reprocessing_MK/poln_calibration/2019-10-07-02:16:00.jones
```

This is difficult for a human to create which is why this file is often created with the [`psrdb`](https://psrdb.readthedocs.io/en/latest/cli.html#observation) command line tool.
For example you can use psrdb to download all observations where the processing failed like so:

```console
psrdb observation download --incomplete --main_project MeerTIME --obs_type fold
```


## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/meerpipe --outdir ./results --pulsar J1410-7404 --utcs 2020-09-05-10:00:47 --utce 2020-09-05-10:00:47 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/meerpipe -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/meerpipe
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/meerpipe releases page](https://github.com/nf-core/meerpipe/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
