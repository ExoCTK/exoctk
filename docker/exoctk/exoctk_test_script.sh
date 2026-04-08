#!/bin/sh

# This script will invoke the various exoctk pages, and test their basic functionality.
# The URL is http://localhost:5000/WHATEVER
# Sub-sites are:
#   /groups_integrations
#       - Enter a target star
#       - Invoke "resolve"
#       - Invoke "calculate"
#   /contam_visibility
#       - Enter a target star
#       - Invoke "resolve"
#       - Invoke "Calculate Visibility"
#       - Enter a target star
#       - Invoke "resolve"
#       - Invoke "Calculate Visibility and Contamination"
#       - Enter a target star
#       - Invoke "resolve"
#       - Set PA to a non-negative-one value
#       - Invoke "Calculate Visibility and Contamination"
#   /limb_darkening
#       - Enter a target star
#       - Invoke "resolve"
#       - Invoke "Calculate Coefficients"
#       - On results page, invoke "Download Coefficients"
#       - On results page, invoke "Download SPAM Coefficients"
#   /fortney
#       - Invoke "submit"
#       - On results page, invoke "Save File"
#   /generic
#       - Set Temperature to 400
#       - Set Gravity to 50
#       - Set planetary radius to 1
#       - Set stellar radius to 1
#       - Invoke "Submit"
#       - On results page, invoke "Save File"
#   /phase_constraint
#       - Enter a target star
#       - Invoke "resolve"
#       - Invoke "Calculate phase constraint"


until cd /exoctk/exoctk/
do
    echo "Waiting for server volume..."
done

# Retrieve the main URL
curl http://localhost:5000
