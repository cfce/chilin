# Django settings for mdseqpos project.

import os, sys
import mdseqpos

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (('Len Taing', 'lentaing@jimmy.harvard.edu'),)

MANAGERS = ADMINS

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/New_York'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# Make this unique, and don't share it with anybody.
SECRET_KEY = ''

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.load_template_source',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
)

ROOT_URLCONF = 'www.urls'

DEPLOY_DIR = mdseqpos.__path__[0]
TEMPLATE_DIRS = (DEPLOY_DIR+'/django')

# THINGS to tailor for your site:
# ASSEMBLY_DIR is where your UCSC genome assemblies are
# BUILD_DICT is a mapping between genome assembly names/short-hands and
# their associated data--note: they must be sub-dirs of ASSEMBLY_DIR
ASSEMBLY_DIR = ''
BUILD_DICT = {
}
