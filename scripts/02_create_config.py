#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all sample names
allSamples = list()
numSamples = 0

with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        numSamples += 1
        line = line.replace(".", "_")
        line = line.replace("-", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # CH_ChMB1_1_BR_Whole_C1_WGMRS_A44752_235MMJLT4_TTGTAAGAGG_L002_R2_001.fastq.gz
                                                #  E1_BR.FCHVC2VDRXY_L1_R1_ITAAGTGGT-CTTAAGCC.fastq.gz pigID_tissue.sequencer_lane_read_X-X.fastq.gz
        # create a shorter sample name
        stemName = sampleAttributes[1] 
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "rawReads" : "/tgen_labs/jfryer/projects/dirty_mice/CH/bulkRNA",
    "rawQC" : "../rawQC/",
    "trimmedReads" : "../trimmedReads/",
    "trimmedQC" : "../trimmedQC/",
    "starAligned" : "../starAligned/",
    "bamstats" : "../bamstats/",

    "Comment_Reference" : "This section specifies the location of the mouse, Ensembl reference genome",
    "Mmusculus.fa" : "/tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/fasta/genome",
    "Mmusculus.gtf" : "/tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/genes/genes",
    "star_ref_index" : "/tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A_STARv2.7.11_150sjdb",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
'''
outfile.write(header.format(allSamples))

# config formatting
counter = 0
with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        base = base.replace("_R1_001", "")
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace(".", "_")
        line = line.replace("-", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # project_uniqueNum_1_tissue_group_XX_XX_sequencer_adapter_lane_read_001.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[1] 
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = stemName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/tgen_labs/jfryer/projects/dirty_mice/CH/bulkRNA/",
        "fq1": "{1}",
        "fq2": "{2}",
        "ID": "{3}",
        "SM": "{4}",
        "PU": "{5}",
        "LB": "{6}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, stemName, stemName, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
