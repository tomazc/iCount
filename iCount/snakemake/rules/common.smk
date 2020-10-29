# Create log folder for cluster run
logdir = os.path.join(os.getcwd(), config["logdir"])
os.makedirs(logdir, exist_ok=True)


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### Completeness Options #####
def all_xlsites(wilcards):
    xlsites_output = list()
    if config["completeness_output"] == "minimum":
        # xlsites_output.extend(expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bed", project=config['project'], barcode=samples.index))
        # xlsites_output.extend(expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed", project=config['project'], barcode=samples.index))
        xlsites_output.extend(expand("{project}/clusters/{barcode}/{barcode}.clusters.bed", project=config['project'],
                                     barcode=samples.index))
    elif config["completeness_output"] == "complete":
        # xlsites_output.extend(expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bed", project=config['project'], barcode=samples.index))
        xlsites_output.extend(
            expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.summary_gene.tsv", project=config['project'],
                   barcode=samples.index))
        xlsites_output.extend(expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_biotype.tab",
                                     project=config['project'], barcode=samples.index))
        xlsites_output.extend(expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.annotated_sites_gene_id.tab",
                                     project=config['project'], barcode=samples.index))
        xlsites_output.extend(
            expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.bedgraph", project=config['project'],
                   barcode=samples.index))
        xlsites_output.extend(
            expand("{project}/xlsites/{barcode}/rnamaps/", project=config['project'], barcode=samples.index))
        # xlsites_output.extend(expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.bed", project=config['project'], barcode=samples.index))
        xlsites_output.extend(
            expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.summary_gene.tsv", project=config['project'],
                   barcode=samples.index))
        xlsites_output.extend(expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.biotype.tab",
                                     project=config['project'], barcode=samples.index))
        xlsites_output.extend(expand("{project}/sig_xlsites/{barcode}/{barcode}.sig_sites.annotated.gene_id.tab",
                                     project=config['project'], barcode=samples.index))
        xlsites_output.extend(expand("{project}/clusters/{barcode}/{barcode}.clusters.bed", project=config['project'],
                                     barcode=samples.index))
    else:
        sys.exit("completeness_output option in configfile needs to be set as minimum or complete ")

    return xlsites_output


def all_group(wilcards):
    group_output = list()
    if config["group_completeness_output"] == "minimum":
        group_output.extend(
            expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed", project=config['project'],
                   group=samples["group"].dropna().unique()))
        group_output.extend(
            expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed", project=config['project'],
                   group=samples["group"].dropna().unique()))
        group_output.extend(
            expand("{project}/groups/{group}/clusters/{group}.group.clusters.bed", project=config['project'],
                   group=samples["group"].dropna().unique()), )
    elif config["group_completeness_output"] == "complete":
        group_output.extend(
            expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bed", project=config['project'],
                   group=samples["group"].dropna().unique()))
        group_output.extend(
            expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_biotype.tab",
                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(
            expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.annotated_group_gene_id.tab",
                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_gene.tsv",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        # group_output.extend(expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_type.tsv", project=config['project'], group=samples["group"].dropna().unique()))
        # group_output.extend(expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.summary_subtype.tsv", project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(
            expand("{project}/groups/{group}/xlsites/{group}.group.unique.xl.bedgraph", project=config['project'],
                   group=samples["group"].dropna().unique()), )
        group_output.extend(expand("{project}/groups/{group}/rnamaps/", project=config['project'],
                                   group=samples["group"].dropna().unique()))

        group_output.extend(
            expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.bed", project=config['project'],
                   group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.biotype.tab",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.annotated.gene_id.tab",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_gene.tsv",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_type.tsv",
                                   project=config['project'], group=samples["group"].dropna().unique()))
        group_output.extend(expand("{project}/groups/{group}/sig_xlsites/{group}.group.sig_sites.summary_subtype.tsv",
                                   project=config['project'], group=samples["group"].dropna().unique()))

        group_output.extend(
            expand("{project}/groups/{group}/clusters/{group}.group.clusters.bed", project=config['project'],
                   group=samples["group"].dropna().unique()), )

    else:
        sys.exit("group_completeness_output option in configfile needs to be set as minimum or complete ")

    return group_output

def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all.
    """
    final_output = []

    if config["bedgraph_UCSC"]:
        final_output.extend(expand("{project}/xlsites/{barcode}/{barcode}.unique.xl.UCSC.bedgraph", project=config['project'], barcode=samples.index))

    if config["create_integrity_test"]:
        final_output.extend(expand("{project}/shasum_files/demultiplex_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/qc_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/genome_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/mapped_reads_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/xlsites_shasum_file.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/group_shasum_file.txt", project=config['project']))

    if config["integrity_test_check"]:
        final_output.extend(expand("{project}/shasum_files/demultiplex_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/qc_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/genome_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/mapped_reads_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/xlsites_shasum_check.txt", project=config['project']))
        final_output.extend(expand("{project}/shasum_files/group_shasum_check.txt", project=config['project']))

    return final_output


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
    return ("{0}/{1}/star_index/SAindex".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_star_index_folder(wildcards):
    return ("{0}/{1}/star_index/".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_segment_path(wildcards):
    return ("{0}/{1}/segment/homo_sapiens_segment.gtf".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_segment_regions(wildcards):
    return ("{0}/{1}/segment/regions.gtf.gz".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

def get_templates_dir(wildcards):
    return ("{0}/{1}/segment/".format(config['genomes_path'], samples.loc[wildcards.barcode, "mapto"]))

