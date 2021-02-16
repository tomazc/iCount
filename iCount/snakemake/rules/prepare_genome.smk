
#==============================================================================#
#                       Download annotation and index genome
#==============================================================================#
# if the genome is not in the config file will fail with a hint to include path to fasta file and annotation.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Check for custom genomes *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Capture iCount available genomes
# --chromosomes 21 MT Include Chromosome from config file!!
# -r 88
args = ["iCount species --source ensembl -r 88"]
# if config["release"]:
#     args.extend(config['release'])
#
# print("RELEASE and download:", args)


species_out = subprocess.Popen(list(args), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
# species_out = subprocess.Popen(["iCount species --source ensembl -r 88"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
stdout,stderr = species_out.communicate()
available_genomes=stdout.decode('utf-8').rstrip()
available_genomes=re.split('available: ',str(available_genomes))[-1].split(',')

all_genomes=samples["mapto"].unique()
custom_genomes=np.setdiff1d(all_genomes, available_genomes)
download_genomes=np.setdiff1d(all_genomes, custom_genomes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* Download ensembl genomes *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tested with homo_sapiens and mus_musculus; custom hg19, hg38, mm10 and mm15.
rule download_genome:
    output:
        genome_fasta="{genomes_path}/{genome}/{genome}.fa.gz",
        genome_index="{genomes_path}/{genome}/{genome}.fa.gz.fai",
        gtf="{genomes_path}/{genome}/{genome}.gtf.gz",
    params:
        release=config['release'],
        chromosomes=config['chromosomes'],
        source=config['source'],
    run:
        GENOME=wildcards.genome
        print ("Adquiring genome: %s \n" % (GENOME))

        if GENOME in download_genomes:
            print ("Downloading iCount available genome:", GENOME)
            print ("Downloading genomes could take some time depending on your conection")
            shell (download_chromosome(output.genome_fasta, params.source, GENOME, params.release, params.chromosomes))
            shell("iCount annotation --annotation {output.gtf} --source {params.source} {wildcards.genome} {params.release}")

        elif GENOME in config['custom_genome'].keys():
            fasta_in = custom_fasta(GENOME)
            gtf_in = custom_gtf(GENOME)

            # Move genome data
            shutil.copy(fasta_in, output.genome_fasta)
            shutil.copy(gtf_in, output.gtf)

            # Create fasta index
            temp = decompress_to_tempfile(fasta_in)
            pysam.faidx(temp)  # pylint: disable=no-member
            shutil.move(temp + '.fai', output.genome_index)

        else:
            print ("Your genome %s in the annotation table %s is not in the iCount available genomes %s \n\n" % (GENOME, config["samples"], available_genomes))
            print ("Please, check misspelled genome or include custom genome %s fasta sequence and annotation GTF file in the config file:" % (GENOME))
            print (yaml.dump(config, default_flow_style=False))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* STAR genome index *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rule indexstar_genome:
    input:
        genome_fasta="{genomes_path}/{genome}/{genome}.fa.gz",
        gtf="{genomes_path}/{genome}/{genome}.gtf.gz",
    threads:
        int(config['index_threads'])
    params:
        overhang=config['overhang'],
    output:
        directory("{genomes_path}/{genome}/star_index/"),
    shell:
        """
        mkdir {output}
        iCount indexstar --overhang {params.overhang} --annotation {input.gtf} \
        --threads {threads} --genome_sasparsed 2 {input.genome_fasta} {output}
        """

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~* iCount genome segment *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rule segment:
    input:
        gtf="{genomes_path}/{genome}/{genome}.gtf.gz",
        genome_fai="{genomes_path}/{genome}/{genome}.fa.gz.fai"
    output:
        segment="{genomes_path}/{genome}/segment/{genome}_segment.gtf",
        landmarks="{genomes_path}/{genome}/segment/landmarks.bed.gz",
    shell:
        """
        iCount segment {input.gtf} {output.segment} {input.genome_fai}
        """
