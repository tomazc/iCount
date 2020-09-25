# Create log folder for cluster run
logdir = os.path.join(os.getcwd(), config["logdir"])
os.makedirs(logdir, exist_ok=True)

# Print project
print("Procesing project:", config['project'], "\n")


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### Helper functions #####

def get_genome(wildcards):
    return ("{0}".format(samples.loc[wildcards.barcode, "mapto"]))

# Custom genomes path from config file
def custom_fasta(wildcards):
    return config['custom_genome'][wildcards]['genome_fasta']

def custom_gtf(wildcards):
    return config['custom_genome'][wildcards]['annotation']


# Function from icount (call it!!)
def decompress_to_tempfile(fname, context='misc'):
    """
    Decompress files ending with .gz to a temporary file and return filename.
    If file does nto end with .gz, juts return fname.
    Parameters
    ----------
    fname : str
        Path to file to open.
    context : str
        Name of temporary subfolder where temporary file is created.
    Returns
    -------
    str
        Path to decompressed file.
    """
    if fname.endswith('.gz'):
        tmp_dir = os.path.join(context)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        suffix = '_{:s}'.format(os.path.basename(fname))
        fout = tempfile.NamedTemporaryFile(suffix=suffix, dir=tmp_dir, delete=False)
        fin = gzip.open(fname, 'r')
        shutil.copyfileobj(fin, fout)
        fin.close()
        fout.close()
        return fout.name

    return fname

def bedgraph_header(wildcards):
    # db=\"{mapto}\" removed
    # return ("{project}_{sample_name}_{protein}_{method}_{mapto}".format(project=config['project'], sample_name=samples.loc[wildcards.barcode, "sample_name"], mapto=samples.loc[wildcards.barcode, "mapto"], method=samples.loc[wildcards.barcode, "method"],	protein=samples.loc[wildcards.barcode, "protein"],	cells_tissue=samples.loc[wildcards.barcode, "cells_tissue"],	condition=samples.loc[wildcards.barcode, "condition"],))
    return ("track type=bedGraph name=\"{project}_{sample_name}_{protein}_{method}_{mapto}_unique.xl.bedgraph.bed\" description=\"{project}_{sample_name}_{protein}_{method}_{mapto}\" "
            "color=\"120,101,172\"  mapped_to=\"{mapto}\" altColor=\"200,120,59\" lib_id=\"{project}\" maxHeightPixels=\"100:50:0\" visibility=\"full\" "
            "tissue=\"{cells_tissue}\" protein=\"{protein}\" species=\"{mapto}\" condition=\"{condition}\" res_type=\"T\" priority=\"20\" \n".format(project=config['project'], sample_name=samples.loc[wildcards.barcode, "sample_name"], mapto=samples.loc[wildcards.barcode, "mapto"], method=samples.loc[wildcards.barcode, "method"],	protein=samples.loc[wildcards.barcode, "protein"],	cells_tissue=samples.loc[wildcards.barcode, "cells_tissue"],	condition=samples.loc[wildcards.barcode, "condition"],))



def get_gtf_path(wildcards):
    return ("{0}/{1}/{1}.gtf.gz".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_star_index_path(wildcards):
    return ("{0}/{1}/star_index/".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_segment_path(wildcards):
    return ("{0}/{1}/segment/{1}_segment.gtf".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_templates_dir(wildcards):
    return ("{0}/{1}/segment/".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

