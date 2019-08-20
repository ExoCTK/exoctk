from distutils.extension import Extension

def get_package_data():
    return {'exoctk': ['data/*', 'data/images/*', 'data/core/*', 'data/core/modelgrid/*', 'data/contam_visibility/*', 'exoctk_app/static/images/*', 'exoctk_app/static/css/*', 'exoctk_app/static/js/*', 'exoctk_app/static/js/vendor/*']}
