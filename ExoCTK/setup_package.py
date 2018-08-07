from distutils.extension import Extension

def get_package_data():
    return {'ExoCTK': ['data/*', 'data/images/*', 'data/core/*', 'data/contam_visibility/*']}
