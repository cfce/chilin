import os
import ConfigParser

#THE following line is cool b/c it allows us to introspect the module!
from . import __path__ as _DEPLOY_DIR

_configParser = ConfigParser.ConfigParser()
_changeset = "Unknown"
if len(_configParser.read(os.path.join(_DEPLOY_DIR[0], "pkg.cfg"))) != 0:
    _changeset = _configParser.get("pkg_info", "hg_changeset")

#NOTE about versions:
#VERSIONS are: preamble (optional) + major version + changeset
#where the preamble and the major version are defined by humans in this file
#and the changeset is automatically updated in the pkg.cfg file by machine
#**THEREFORE: you should only need to modify _preamble and _major_version.
_preamble = "mdseqpos (official trunk)"
_major_version = "2.01"
__version__ = "%s: Version %s\nChangeset: %s" % (_preamble, _major_version, _changeset)
