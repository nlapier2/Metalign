import logging
import re
import sys
from subprocess import Popen, PIPE


def check_minimap2_bin(minimap2_bin='minimap2', version=r'^2\.\S+'):
    """Check minimap2 binary exists and is the correct version

    Running
    $ minimap2 --version
    should produce a single line with just the version of minimap2, e.g "2.17-r941"

    Log issue and exit immediately if correct version of minimap2 cannot be found.

    Args:
        minimap2_bin (str): Minimap2 binary path
        version (str): Minimap2 version regex pattern

    Returns:
        Valid minimap2 path string

    Raises:
        ValueError: when unknown version of minimap2 is specified
        FileNotFoundError: when minimap2 cannot be found at the specified path
    """
    try:
        p = Popen([minimap2_bin, '--version'], stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()
        minimap2_version = stdout.decode().strip()
        m = re.match(version, minimap2_version)
        if m:
            logging.info(f'Found minimap2 version {minimap2_version}')
            return minimap2_bin
        else:
            error_msg = f'Found unsupported version ({minimap2_version}) of minimap2'
            logging.error(error_msg)
            sys.exit(1)
    except FileNotFoundError:
        error_msg = f'Could not find "minimap2" at "{minimap2_bin}". Please install minimap2. e.g. with Conda: $ ' \
                    f'conda install -c bioconda minimap2'
        logging.error(error_msg)
        sys.exit(1)


def check_kmc_bin(kmc_bin='kmc', version=r'\b3\.\S+\b'):
    """Check KMC binary exists and is the correct version (3.X.X)

    Also checks that `kmc_tools` and `kmc_dump` exist and are the correct version.

    The first line of running `kmc` should contain something similar to:

    "K-Mer Counter (KMC) ver. 3.1.1 (2019-05-19)"

    Default `version` pattern is looking for the "3." and non-whitespace characters.

    Log issue and exit immediately if correct version of kmc cannot be found.

    Args:
        kmc_bin (str): KMC3 binary path
        version (str): KMC3 version regex pattern

    Returns:
        Valid kmc path string

    Raises:
        ValueError: when unknown version of kmc is specified
        FileNotFoundError: when kmc cannot be found at the specified path
    """
    try:
        for exe in [kmc_bin, kmc_bin + '_dump', kmc_bin + '_tools']:
            p = Popen([exe, ], stdout=PIPE, stderr=PIPE)
            stdout, _ = p.communicate()
            ver = stdout.decode().strip().splitlines()[0]
            m = re.search(version, ver)
            if m:
                logging.info(f'Found {exe} version {ver}')
            else:
                error_msg = f'Found unsupported version ({ver}) of kmc at "{exe}". ' \
                            f'Please install a supported version (3.X.X)'
                logging.error(error_msg)
                sys.exit(1)
        return kmc_bin
    except ValueError:
        error_msg = f'Unknown version of kmc installed at {kmc_bin}. ' \
                    f'Please installed a supported version (3.X.X)'
        logging.error(error_msg)
        sys.exit(1)
    except FileNotFoundError:
        error_msg = f'Could not find "kmc" at "{kmc_bin}". Please install kmc. e.g. with Conda: $ ' \
                    f'conda install -c bioconda kmc'
        logging.error(error_msg)
        sys.exit(1)
