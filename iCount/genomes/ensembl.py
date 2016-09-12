"""
Functions to querry and download data from ENSEMBL ftp site
"""

import os
import re
import ftplib
import gzip
import shutil
import tempfile
import subprocess

BASE_URL = 'ftp.ensembl.org'

MIN_RELEASE_SUPPORTED = 59
MAX_RELEASE_SUPPORTED = 84


def _docstring_parameter(*sub):
    """
    To be used as decorator for passing arguments for docstring formatting.
    """
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj
    return dec


def get_ftp_instance():
    """
    Get ftplib.FTP object that is connected to BASE_URL

    :return: FTP object connected to BASE_URL
    :rtype: ftplib.FTP
    """
    ftp = ftplib.FTP(BASE_URL)
    ftp.login()
    return ftp


@_docstring_parameter(MIN_RELEASE_SUPPORTED, MAX_RELEASE_SUPPORTED)
def get_release_list():
    """
    Get list of available ENSEMBL releases.

    Only allows ENSEMBL releases from {0} - {1}.

    Returns
    -------
    list
        List of available releases

    """
    ftp = get_ftp_instance()
    # set current working directory
    ftp.cwd('pub')

    out = [item.strip('release-') for item in ftp.nlst() if
           re.match(r'release-\d+', item) and
           int(item.strip('release-')) >= MIN_RELEASE_SUPPORTED and
           int(item.strip('release-')) <= MAX_RELEASE_SUPPORTED]

    ftp.quit()
    return sorted(out, reverse=True)


def get_species_list(release):
    """
    Get list of species for given release.

    Parameters
    ----------
    release : str
        The release number of type str or int.

    Returns
    -------
    list
        List of species.

    """
    ftp = get_ftp_instance()
    ftp.cwd('pub/' + 'release-' + str(release) + '/fasta/')
    species_list = [item for item in ftp.nlst()]

    ftp.quit()
    if species_list:
        return sorted(species_list)


def download_annotation(release, species, target_dir=None, target_fname=None):
    """
    Download annotation in GTF file fromat.

    Parameters
    ----------
    release : int
        Release number.
    species : str
        Species latin name.
    target_dir : str
        Download to this directory (if not given, current working directory).
    target_fname : str
        Annotation filename (must have .gz file extension). If not given,
        species.release.gtf.gz is used.

    Returns
    -------
    str
        Downloaded annotation filename.

    TODO: target_dir & target_fname into one parameter? But the default name is
    a quite useful feature...
    """

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))

    ftp = get_ftp_instance()
    server_dir = '/pub/release-{}/gtf/{}/'.format(release, species)
    ftp.cwd(server_dir)
    server_files = ftp.nlst()

    # Get the desired annotation file form list of all server files:
    annotation_file = ''
    regex = r'{}\.[\w\d]+\.\d+\.gtf\.gz'.format(species.capitalize())
    for file_ in server_files:
        if re.match(regex, file_):
            annotation_file = file_

    if not annotation_file:
        return None

    # Process filename
    if not target_fname:
        target_fname = '{}.{}.gtf.gz'.format(species, release)

    saveas_fname = os.path.join(target_dir, target_fname)

    # Download to file on disk
    with open(saveas_fname, 'wb') as fhandle:
        ftp.retrbinary('RETR ' + annotation_file, fhandle.write)

    ftp.quit()
    return saveas_fname


def chrom_length(fasta_in, txt_out=None):
    """
    Compute chromosome lengths to a file

    Output file has one line per each chromosome::

        chr1    249250621
        chr2    243199373
        ...

    :param string fasta_in: path to genome fasta file (can be *.gz file)
    :param string txt_out: output *.txt file
    :return: absoulute path to output file
    :rtype: str
    """
    if fasta_in.endswith('.gz'):
        f2 = tempfile.NamedTemporaryFile(delete=False)
        with gzip.open(fasta_in, 'rb') as f1:
            shutil.copyfileobj(f1, f2)
        f2.close()
        temp = f2.name
    else:
        temp = fasta_in

    command = ['samtools', 'faidx', temp]
    subprocess.check_call(command)
    # This command makes fai file:
    fai_file = temp + '.fai'

    if not txt_out:
        txt_out = fasta_in + '.chrom_length.txt'

    with open(fai_file, 'r') as f1, open(txt_out, 'wt') as f2:
        for line in f1:
            chrom, length = line.strip().split()[:2]
            f2.write('{}\t{}\n'.format(chrom, length))

    # Clean up:
    os.remove(fai_file)

    return os.path.abspath(txt_out)


def download_sequence(release, species, target_dir=None, target_fname=None,
                      tempdir=None, chromosomes=[]):
    """
    Downloads genome file in FASTA fromat.

    Several steps are performed:

        * querry for list off all FASTA files for given release and specias
        * filter this list to get only whole chromosome files
        * if chromosomes paramter is given, take only specified chromosomes
        * sort list of these files to have the correct order
        * download each file and write it to target_fname

    Parameters
    ----------
    release : int
        Release number.
    species : str
        Species latin name.
    target_dir : str
        Download to this directory (if not given, current working directory).
    target_fname : str
        Annotation filename (must have .gz file extension). If not given,
        species.release.gtf.gz is used.
    tempdir : str
        Temporary directory with intermediate results.
    chromosomes : list_str
        If given, don't download the whole genome, but juts the given
        cromosomes. Chromosomes can be given as strings aor integers

    Returns
    -------
    str
        Downloaded genome/sequnce filename.

    TODO: target_dir & target_fname into one parameter? But the default name is
    a quite useful feature...
    """

    if not target_dir:
        target_dir = os.getcwd()
    if not os.path.isdir(target_dir):
        raise ValueError('Directory "{}" does not exist'.format(target_dir))
    if not target_fname:
        target_fname = '{}.{}.fa.gz'.format(species, release)

    ftp = get_ftp_instance()
    server_dir = '/pub/release-{}/fasta/{}/dna'.format(release, species)
    ftp.cwd(server_dir)
    all_fasta_files = ftp.nlst()

    filtered_files = []

    # Recognize all "whole-chromosme" files
    regex = r'{}\.[\w\d]+\.(\d+\.)*dna\.chromosome\.' \
            r'[\dXYMT]+\.fa\.gz'.format(species.capitalize())
    for fasta_file in all_fasta_files:
        if re.match(regex, fasta_file):
            filtered_files.append(fasta_file)

    def sorting_func(name):
        """Helper function for sorting files named by chromosome"""
        key = name.split('.')[-3]
        if key == 'MT':
            return 'ZZ'  # this makes mitohondrial "chromosome" the last one
        return key
    # Sort the filtered_files in correct order of chromosomes:
    filtered_files.sort(key=sorting_func)

    # If parameter chromosomes is given, take only a subset of all chromosomes
    if chromosomes:
        subset_list = []
        for file_ in filtered_files:
            for chromosome in chromosomes:
                regex = r'.*chromosome\.{}\.fa\.gz'.format(str(chromosome))
                if re.match(regex, file_):
                    subset_list.append(file_)
                else:
                    pass
        filtered_files = subset_list
        target_fname = target_fname.rstrip('.fa.gz') + '.chr' + \
            '_'.join(map(str, chromosomes)) + '.fa.gz'

    tempdir = tempfile.mkdtemp(dir=tempdir)
    target_path = os.path.join(target_dir, target_fname)

    for fname in filtered_files:
        # Download all files to tempdir:
        temp_file = os.path.join(tempdir, fname)
        with open(temp_file, 'wb') as fhandle:
            ftp.retrbinary('RETR ' + fname, fhandle.write)

        # Write the content into final file (target_fname)
        with gzip.open(temp_file, 'rb') as f_in, \
             gzip.open(target_path, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Clean up: delete temporary files:
    ftp.quit()
    shutil.rmtree(tempdir)

    # Conpute chromosome lengths:
    _ = chrom_length(target_path)

    return target_path


# #############################################################

if __name__ == '__main__':

    releases = get_release_list()
    release84 = releases[-2]
    species = get_species_list(releases[-2])

    # Choose random species
    species = species[42]

    # download_annotation(release84, species, '/home/jure/Desktop')

    # Executing this command might take a looong time!
    # download_genome(release84, species, '/home/jure/Desktop')
