"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'iCount'
__summary__ = 'Computational pipeline for analysis of iCLIP data'
__url__ = 'https://github.com/tomazc/iCount'

# Semantic versioning is used. For more information see:
# https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
__version__ = '2.0.0-dev'

__author__ = 'University of Ljubljana, Bioinformatics Laboratory'
__email__ = 'tomazc@gmail.com'

__license__ = 'MIT'
__copyright__ = '2016, ' + __author__

__all__ = (
    '__title__', '__summary__', '__url__', '__version__', '__author__',
    '__email__', '__license__', '__copyright__',
)
