#!/bin/bash

wd=${PWD}
script_dir=$(dirname `readlink -f $0`)
proactive_workflow_xml="${script_dir}/proactive_eskrim_workflow.xml"
proactive_client_exec="/opt/proactive/default/bin/proactive-client"
proactive_default_credentials="${HOME}/.proactive/credentials.txt"
eskrim_path="/usr/local/bin/eskrim.py"

usage="
USAGE: $0 [OPTIONS]

OPTIONS

  -k PROACTIVE_CREDENTIALS             Optional, path to user's credentials file (for proactive submission)
                                       Default: ${proactive_default_credentials}

  -e ESKRIM_PATH                       Optional, path to the ESKRIM script
                                       Default: ${eskrim_path}

  -s SAMPLES_DIR                       MANDATORY, path to folder containing sample directories with fastq files
                                       e.g. /projects/prj_name/sample

  -t TEMP_OUTPUT_DIR                   MANDATORY, directory where temporary results will be stored
                                       e.g. /projects/prj_name/temp_analysis

  -o OUTPUT_FILE                       MANDATORY, output file where final results will be written
                                       e.g. /projects/prj_name/analysis/eskrim_output.tsv
                             
  -1 FORWARD_READS_SUFFIX              Optional, suffix of fastq files corresponding to forward reads
                                       Default: None

  -r READS_LENGTH                       Optional, discard reads shorter than READS_LENGTH bases and trim those exceeding this length 
                                       Default: 80

  -k KMERS_LENGTH                      Optional, length of kmers to count
                                       Default: 21

  -n NUM_READS                         Optional, NUM_READS to draw randomly from fastq files
                                       Default: 10000000

  -f FORCE_OVERWRITE                   Optional, force overwrite results of previous run
                                       Default: False

  -h                                   Print this help message and exit
"

[[ -z $1 ]] && echo "$usage" && exit 0

####################################################################
########### parsing arguments with getopts (bash built-in) #########
####################################################################
while getopts ':k:e:s:t:o:1:r:k:n:f:h' flag; do
        # if argument is required and not given then flag is set to ':'
        #~ [[ ${OPTARG%-[cwlsptMPTNLO]} != $OPTARG ]] && OPTARG=$flag && flag=:
        case $flag in
                k)  PROACTIVE_CREDENTIALS="$OPTARG" ;;
                e)  ESKRIM_PATH="$OPTARG" ;;
                s)  SAMPLES_DIR="$OPTARG" ;;          # Mandatory
                t)  TEMP_OUTPUT_DIR="$OPTARG" ;;      # Mandatory
                o)  OUTPUT_FILE="$OPTARG" ;;          # Mandatory
                1)  FORWARD_READS_SUFFIX="$OPTARG" ;;
                r)  READS_LENGTH="$OPTARG" ;;
                k)  KMERS_LENGTH="$OPTARG" ;;
                n)  NUM_READS="$OPTARG" ;;
                f)  FORCE_OVERWRITE="$OPTARG" ;;
                h)  echo "$usage" && exit 0 ;;
                \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
                :)  echo "argument required for option -$OPTARG" >&2; exit 1 ;;
                *) echo "Unimplemented option: -$OPTARG" >&2; exit 1 ;;
        esac
done

shift $(( OPTIND - 1 ));

##### ARGUMENTS CHECKING #####

### check mandatory options
fail=
[[ -z $SAMPLES_DIR ]] && fail="${fail} -s"
[[ -z $TEMP_OUTPUT_DIR ]] && fail="${fail} -t"
[[ -z $OUTPUT_FILE ]] && fail="${fail} -o"
[[ -n $fail ]] && echo "Error, missing mandatory option/argument(s): ${fail}" >&2 && exit 1

SAMPLES_DIR=$(realpath ${SAMPLES_DIR})
TEMP_OUTPUT_DIR=$(realpath ${TEMP_OUTPUT_DIR})
OUTPUT_FILE=$(realpath ${OUTPUT_FILE})
prj_name=$(basename $(dirname ${SAMPLES_DIR}))

### check and set optional options
[[ -z $PROACTIVE_CREDENTIALS ]] && PROACTIVE_CREDENTIALS=${proactive_default_credentials}
[[ -z $ESKRIM_PATH ]] && ESKRIM_PATH=${eskrim_path}
[[ -z $FORWARD_READS_SUFFIX ]] && FORWARD_READS_SUFFIX=''
[[ -z $READS_LENGTH ]] && READS_LENGTH='80'
[[ -z $KMERS_LENGTH ]] && KMERS_LENGTH='21'
[[ -z $NUM_READS ]] && NUM_READS='10000000'
[[ -z $FORCE_OVERWRITE ]] && FORCE_OVERWRITE='False'

proactive_workflow_params="\
ESKRIM_PATH            ${ESKRIM_PATH}
SAMPLES_DIR            ${SAMPLES_DIR}
TEMP_OUTPUT_DIR        ${TEMP_OUTPUT_DIR}
OUTPUT_FILE            ${OUTPUT_FILE}
FORWARD_READS_SUFFIX   ${FORWARD_READS_SUFFIX}
READS_LENGTH           ${READS_LENGTH}
KMERS_LENGTH           ${KMERS_LENGTH}
NUM_READS              ${NUM_READS}
FORCE_OVERWRITE        ${FORCE_OVERWRITE}"

proactive_wokflow_params_short="{`echo -ne "$proactive_workflow_params" | perl -ne 's/(\S+)/\"$1\"/g; s/ +/:/g; s/\n/, /g; print'`}"
proactive_wokflow_params_short=$(echo ${proactive_wokflow_params_short} | sed 's/":,/":"",/')
tmp_proactive_workflow_xml=$(mktemp --suffix '.xml')
sed -e "s/eskrim prj_name/eskrim ${prj_name}/" $proactive_workflow_xml > $tmp_proactive_workflow_xml

echo
echo "submitting fastp job with following parameters:"
echo
echo "$proactive_workflow_params"

$proactive_client_exec -X -k -u http://skoody/rest -c ${PROACTIVE_CREDENTIALS} -s ${tmp_proactive_workflow_xml} "${proactive_wokflow_params_short}"

exit 0
