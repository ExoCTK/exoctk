#!/bin/sh

until cd /exoctk/exoctk/
do
    echo "Waiting for server volume..."
done

# run gunicorn
gunicorn -b 0.0.0.0:5000 --log-level info --access-logfile - --error-logfile - --capture-output 'exoctk.exoctk_app.app_exoctk:app_exoctk'

tail -f /dev/null
