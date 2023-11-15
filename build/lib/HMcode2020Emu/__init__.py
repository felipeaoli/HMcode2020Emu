import numpy as np
import copy
import os
from .matter_powerspectrum import *
import cosmopower

base_dir = os.path.join(os.path.dirname(__file__))


def download_data(download_dir):
    """
    Download the data needed for the emulators to the specified directory.

    Parameters
    ----------
    download_dir : str
        the data will be downloaded to this directory
    """
    from six.moves import urllib
    import shutil
    import gzip
    import glob

    #url = 'https://drive.google.com/file/d/1ONzDCLQQ3N6LH_IfXASnP3glH9XlE0yX/view?usp=sharing'
    filenames = 'models.zip'
    out_filenames = 'models'
    file_path = os.path.join(download_dir, filenames)
    final_path = os.path.join(download_dir, out_filenames)
    import gdown

    url = 'https://drive.google.com/uc?id=1ONzDCLQQ3N6LH_IfXASnP3glH9XlE0yX'


    # do not re-download
    if not os.path.exists(final_path):

        print("\n As it is the first instance of the emulator, "
                    "we need to download some data, it can take a few "
                    "seconds...\n")

        print("Downloading %s...\n" % out_filenames)

        #file_path, _ = urllib.request.urlretrieve(url=url,
        #                                            filename=file_path,
        #                                            reporthook=None)
        gdown.download(url, file_path, quiet=False)

        print("Download finished. Extracting files.")

        # unzip the file
        shutil.unpack_archive(
            filename=file_path, extract_dir=download_dir)
        os.remove(file_path)
        print("Done.\n")


download_data(base_dir)
