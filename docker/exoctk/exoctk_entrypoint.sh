#!/bin/sh

until cd /exoctk/exoctk/exoctk/exoctk_app
do
    echo "Waiting for server volume..."
done

# run flask
# flask --app app_exoctk run

python /exoctk/exoctk/exoctk/exoctk_app/exoctk_app.py
