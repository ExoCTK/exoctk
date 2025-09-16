#!/bin/sh

until cd /exoctk/exoctk
do
    echo "Waiting for server volume..."
done

# run flask
# flask --app app_exoctk run

# python /exoctk/exoctk/exoctk/exoctk_app/app_exoctk.py

# run gunicorn
gunicorn -w 4 -b 0.0.0.0 'exoctk.exoctk_app.app_exoctk:app_exoctk'
