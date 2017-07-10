from distutils.extension import Extension

def get_package_data():
    return {'ExoCTK': ['data/*', 'data/filters/*', 'data/tot/*', 'data/pal/*', 'data/bar/*', 'data/ldc/*', 'data/core/*', 'data/tor/*']}
