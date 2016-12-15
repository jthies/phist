#!/bin/bash
export RM_SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pbsdsh -u ${RM_SCRIPT_PATH}/do-rm-cache-ckpts.sh
