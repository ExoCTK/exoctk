from distutils.extension import Extension

def get_package_data():
    return {'ExoCTK': ['data/*', 'data/filters/*']}
