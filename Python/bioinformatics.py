def blast_import(file_path):
    '''
    imports individual blast text file and adds column names
    file_path - path to text file
    '''
    dataframe = pd.read_table(file_path, sep = "\t", header = None)
    blast_colnames = ["query seqid", "subject seqid", "identity match %", "length", "mismatches", "gap open",
                  "qstart", "qend", "alignment start", "alignment end",
                  "evalue", "bitscore"]
    dataframe.columns = blast_colnames
    return(dataframe)

def blast_columns(dataframe):
    '''
    adds blast column names to blast dataframe from text file
    must import dataframe first from text - more limited than blast_import
    '''
    if isinstance(dataframe,pd.DataFrame) == True:
        blast_colnames = ["query seqid", "subject seqid", "identity match %", "length", "mismatches", "gap open",
                      "qstart", "qend", "alignment start", "alignment end",
                      "evalue", "bitscore"]
        dataframe.columns = blast_colnames
        return(dataframe.head())
    else:
        dataframe = None
        print("Provide a pandas dataframe")

def fasta_concat(directory, file_end, out_dir_and_name):
    '''
    creates one fasta file from many individual fasta files
    useful for creating a custom blast db

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
