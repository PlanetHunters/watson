__version__ = "1.1.0"

# Patching all errors due to SSL certificates + TLS version fix
import requests
import warnings
import urllib3
warnings.filterwarnings("ignore", category=urllib3.exceptions.InsecureRequestWarning)

import ssl

# MAST only supports TLS 1.2. OpenSSL 3.6+ defaults to TLS 1.3 which causes
# [SSL] record layer failure. Limit all SSL contexts to TLS 1.2 maximum.
def _create_tls12_unverified_context():
    ctx = ssl._create_unverified_context()
    ctx.minimum_version = ssl.TLSVersion.TLSv1_2
    ctx.maximum_version = ssl.TLSVersion.TLSv1_2
    ctx.set_ciphers('DEFAULT@SECLEVEL=1')
    return ctx

ssl._create_default_https_context = _create_tls12_unverified_context

# Patch urllib3's context creator to also limit TLS version
_original_create_urllib3_context = urllib3.util.ssl_.create_urllib3_context

def _patched_create_urllib3_context(ssl_version=None, cert_reqs=None, options=None, ciphers=None,
                                     ssl_minimum_version=None, ssl_maximum_version=None, verify_flags=None):
    if ssl_maximum_version is None:
        ssl_maximum_version = ssl.TLSVersion.TLSv1_2
    if ssl_minimum_version is None:
        ssl_minimum_version = ssl.TLSVersion.TLSv1_2
    return _original_create_urllib3_context(
        ssl_version=ssl_version, cert_reqs=cert_reqs, options=options,
        ciphers=ciphers, ssl_minimum_version=ssl_minimum_version,
        ssl_maximum_version=ssl_maximum_version, verify_flags=verify_flags
    )

urllib3.util.ssl_.create_urllib3_context = _patched_create_urllib3_context

# Original request method
_original_request = requests.Session.request
# Overwrite to enforce verify=False by default
def no_verify_request(self, *args, **kwargs):
    kwargs.setdefault("verify", False)  # solo cambia si no se pasa ya 'verify'
    return _original_request(self, *args, **kwargs)
# Apply patch
requests.Session.request = no_verify_request

from watson.watson.watson import Watson