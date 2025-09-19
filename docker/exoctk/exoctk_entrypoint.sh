#!/bin/sh

until cd /exoctk/exoctk/
do
    echo "Waiting for server volume..."
done

# run flask
# flask --app app_exoctk run

# run gunicorn
gunicorn -b 0.0.0.0:5000 --log-level debug 'exoctk.exoctk_app.app_exoctk:app_exoctk'
