## BLAST ##

def blast_import(file_path, delim):
    '''
    imports individual blast text file and adds column names
    inputs:
        file_path - path to text file
        delim - type of delimiter e.g. "\t", ","
    '''
    import pandas as pd

    dataframe = pd.read_table(file_path, sep = delim, header = None)
    blast_colnames = ["query_seqid", "subject_seqid", "identity_match", "length", "mismatches", "gap_open",
                  "qstart", "qend", "alignment_start", "alignment_end",
                  "evalue", "bitscore"]
    dataframe.columns = blast_colnames
    return(dataframe)

def blast_concat(file_end, out_name):
    '''
    imports multiple blast files, such as individual files from blasting MAGs
    returns all the data in one pandas dataframe with a column for the file name
    assumes files are in the same folder as the script

    inputs:
        file_end - common file ending to search for such as "*_blast_cphA.txt"
            ideally file names should a similar structure to this: seqname_searchtype_genequery.txt
        out_name - output name, will save to current working directory unless another is specified
    example:
        blast_concat("*_blast_cphA.txt","all_blast.csv")

    '''
    import glob
    import csv

    alldata_list = []

    filenames = glob.glob(file_end)

    for name in filenames:
        with open(name) as read_file:
            reader = csv.reader(read_file, delimiter="\t")
            file_list = list(reader)
            temp_end = file_end.replace("*","") # string with file ending to remove
            seqid = name.replace(temp_end,"")
        for l in file_list:
            l.append(seqid)
            alldata_list.append(l)
    with open(out_name, "w", newline="") as write_file:
        writer = csv.writer(write_file)
        writer.writerows(alldata_list)

def blast_columns(dataframe):
    '''
    adds blast info column names to blast dataframe from text file
    must import dataframe first from text - more limited than blast_import
    '''
    import pandas as pd

    if isinstance(dataframe,pd.DataFrame) == True:
        blast_colnames = ["query_seqid", "subject_seqid", "identity_match", "length", "mismatches", "gap_open",
                  "qstart", "qend", "alignment_start", "alignment_end",
                  "evalue", "bitscore"]
        dataframe.columns = blast_colnames
        return(dataframe.head())
    else:
        dataframe = None
        print("Provide a pandas dataframe")

def fasta_concat(directory, file_end, out_dir_and_name):
    '''
    creates one fasta file from many individual fasta files
    useful for creating a custom blast db or combining multiple fa files

    inputs:
        directory - location where individual fasta files are saved
        file_end - common file name ending
        out_dir_and_name - output directory and file name

    example:
        fasta_concat("./Sequences", "/*cphA.fasta", "./Database files/cphA.fasta")
    '''
    import glob
    from Bio import SeqIO

    filenames = glob.glob(directory + file_end)
    with open(out_dir_and_name, "w") as w_file:
        for file in filenames:
            record = SeqIO.parse(file, "fasta")
            SeqIO.write(record, w_file, "fasta")

## CHECKM ##
def checkm_qa_parse(file_path,out_path):
    '''
    creates a json file from a checkM qa file (https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics#qa)
    the qa file is a tsv, which can be hard to work with directly
    makes a nested dictionary, where the format is {bin id: {qa statistic: value, ...}}

    inputs:
        file_path = where the tsv file is
        out_path = location and name of output file, make sure this is a .json file

    '''
    import json

    path = file_path
    qa = open(path).read()
    # cleaning - parsing the file into a list of lines
    qa_lines = qa.split("\n")
    # removing last list entry if its blank
    if qa_lines[-1] == "":
        qa_lines.pop()
    else:
        pass
    # save to dict
    qa_dict = {}
    for line in qa_lines:
        line_temp = line.split("\t")
        data = json.loads(line_temp[1].replace("'",'"'))
        qa_dict[line_temp[0]] = data
    # export to json
    with open("checkm_qa.json", "w", encoding='utf-8') as f:
        json.dump(qa_dict, f, ensure_ascii=False)

## FILE COPY/MOVE ##
def file_copy(in_dir, out_dir, list_name):
    """
    copy files from individual folders to one folder
    inputs:
        in_dir = path to find files, if file is saved in prokka dir, use "./*/*.tsv"
        out_dir = where files should go, this folder has to exist before running the code!!
        list_name = a name to keep a list of files that were copied
    """
    import glob
    import json
    from shutil import copy

    files = glob.glob(in_dir)

    with open(list_name, "w") as outfile:
        json.dump(files, outfile) # write a json output

    for i in files:
        copy(i,out_dir)

def file_move(in_dir, out_dir, list_name):
    """
    move files from individual folders to one folder
    inputs:
        in = path to find files, if file is saved in prokka dir, use "./*/*.tsv"
        out_dir = where tsv files should go, this folder has to exist before running the code!!
        list_name = a name to keep a list of files that were moved
    """
    import glob
    import json
    from shutil import move

    files = glob.glob(in_dir)

    with open(list_name, "w") as outfile:
        json.dump(files, outfile) # write a json output

    for i in files:
        move(i,out_dir)

## PROKKA ##
def prokka_search(directory, search_term):
    '''
    search tsvs generated by prokka for a particular search term such as a gene of interest
    inputs:
        directory = where the tsv files are located
        search_term = what you're looking for; example - "cphA"
    outputs:
        result_df = pandas dataframe of prokka info and bin
    '''
    import glob
    import pandas as pd

    match_list = []
    line_list = []
    replace_string = directory.replace("/*.tsv","") + "\\"
    filenames = glob.glob(directory)

    for name in filenames:
        with open(name) as read_file: # open each file
            all_lines = read_file.readlines() # read the lines in
            for line in all_lines:
                if search_term in line: # find search term in line
                    match_list.append(name)
                    line_list.append(line)
                else:
                    pass

    for i in range(0,len(match_list)):
        match_list[i] = match_list[i].replace(replace_string,"")
        match_list[i] = match_list[i].replace(".tsv","")
        line_list[i] = line_list[i].replace("\n","")
    result_df = pd.DataFrame([x.split("\t") for x in line_list], columns=["locus_tag","ftype","length_bp","gene","EC_number","COG","product"])
    result_df["bin"] = match_list

    return result_df

def prokka_gff_parse(directory_gff,target_ec):
    '''
    Parsing through prokka-generated gff files to get contig information for annotations of interest
    inputs:
        directory_gff = where gff files are located, should also include /*.gff to be safe
        target_ec = target enzyme commission number
        output = file name
    outputs:
        target_contigs_df = df of bins/contigs and annotations
        need to set to a variable name
    '''
    import glob
    import pandas as pd

    filenames = glob.glob(directory_gff) # list of filenames to iterate through
    target_contigs = [] # empty list where ALL results will go

    column_names = ["bin:contig","start","end","locus_tag","ec_num"]
    keep_index = [0,3,4,8,9] # which target_contigs_columns to remove
    # this is dependent on the prokka version and gff format
    # I am trying to keep bin/contig ID, begin and end of sequence, prokka-assigned ID, EC number

    for name in filenames: # outer loop layer - open each gff file and perform these functions
        # append every thing to a common list/df
        clean_lines = [] # empty list where contig info lines will go
        target_contig_list = [] # empty list where target contig lines will go each loop

        with open(name) as read_file: # open each file, will loop with outer loop
            all_lines = read_file.readlines() # read lines in - splits by new line
            clean_lines = [i for i in all_lines if not ("##" in i or ">" in i)] # removes conttig sequence lines
            clean_lines = [i.strip() for i in clean_lines] # removes \n character

            for i in clean_lines: # inner loop, looking for specific ec number in each line
                if "eC_number="+target_ec in i:
                    target_contig_list.append(i) # appending entry to list

            # cleaning lines - had to split for this to work on quest for some unknown reason
            target_contig_list = [item.replace("\t",";") for item in target_contig_list] # cleaning line
            target_contig_list = [item.split(";") for item in target_contig_list] # cleaning line

            target_contigs.extend(target_contig_list) # append all results to a list

            # # removing substrings from remaining data since this will go into a dataframe with column names
            target_contigs = [[x.replace("ID=","") for x in item] for item in target_contigs]
            target_contigs = [[x.replace("eC_number=","") for x in item] for item in target_contigs]
            target_contigs = [[x.replace("product=","") for x in item] for item in target_contigs]

    # convert target_contigs list into dataframe for downstream use
    target_contigs_df = pd.DataFrame(target_contigs)
    target_contigs_df = target_contigs_df[keep_index] # only keeping relevant columns listed above in keep_index
    target_contigs_df.columns = column_names

    return target_contigs_df

## QIIME
def qiime_manifest(directory, expression, output):
    '''
    generate a qiime manifest file from a directory of 16s data
    inputs:
        directory = where the 16s files are located, should be in fastq.gz format
        expression = regex expression of unique sample ids
        output = name of manifest file
    outputs:
        manifest csv in same directory as script
    '''
    import glob
    import pandas as pd
    import re

    directory = directory + "*.fastq.gz"
    filenames = glob.glob(directory)

    data = [] # empty list to fill, eventually make into DataFrame
    sample_id = []
    filepath =[]
    direction = []

    forward = "_R1_"
    reverse = "_R2_"

    for i in filenames: # sample ids
        sample_id.append(re.search(expression,i).group())

    for i in filenames: # read direction
        if forward in i:
            direction.append("forward")
        elif reverse in i:
            direction.append("reverse")

    data = [sample_id,filenames,direction] # combining lists
    df = pd.DataFrame(data).transpose() # need to transpose list of lists
    df.columns = ["sample-id","absolute-filepath","direction"]
    df.to_csv(output, index=False) #export
