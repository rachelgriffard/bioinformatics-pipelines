
# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions
CONDA=/kuhpc/work/biostat/r816g589/conda/envs/

#THIS MUST BE AT THE END OF THE FILE FOR SDKMAN TO WORK!!!
export SDKMAN_DIR="/home/r816g589/.sdkman"
[[ -s "/home/r816g589/.sdkman/bin/sdkman-init.sh" ]] && source "/home/r816g589/.sdkman/bin/sdkman-init.sh"
