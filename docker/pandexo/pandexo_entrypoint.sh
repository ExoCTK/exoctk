#!/bin/sh

until cd /pandexo
do
    echo "Waiting for server volume..."
done

# Start the first process in the background
start_pandexo --workers=5 &

# Wait for any background process to exit
wait -n

# Exit with the status of the process that exited first
exit $?
