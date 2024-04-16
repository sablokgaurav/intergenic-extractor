def generateintergenic(inputgff, fasta, outfile):
    """
    # Author: Gaurav Sablok
    # Universitat Potsdam
    # Date: 2024-4-16
    a integenic regions stichter for the genome annotations and can 
    annotated and extract all the integenic regions present by the alignment
    of the protein to the genome regions.
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".coding.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    read_transcripts = [i.strip() for i in open(fasta, "r").readlines()]
    fasta_transcript_dict = {}
    for i in read_transcripts:
        if i.startswith(">"):
            path = i.strip()
            if i not in fasta_transcript_dict:
                fasta_transcript_dict[i] = ""
                continue
        fasta_transcript_dict[path] += i.strip()
    fasta_sequences = list(fasta_transcript_dict.values())
    fasta_names = [i.replace(">", "") for i in list(fasta_transcript_dict.keys())]
    readiterator = [i.strip().split() for i in open(inputgff + ".coding.gff").readlines() if i.strip().split()[2] == "CDS"]
    ids = set([read[i][0] for i in range(len(readiterator))])
    intergenic = []
    for i in range(len(read)):
        if read[i][0] in ids:
            intergenic.append([read[i][0], read[i][3], read[i][4]])
        else:
            pass
    extractintergenic = []
    for i in range(len(intergenic)-1):
        for j in range(len(fasta_names)):
            for k in range(len(fasta_sequences)):
                if intergenic[i][0] == fasta_names[j]:
                    extractintergenic.append([intergenic[i][0], fasta_sequences[k][int(intergenic[i][2]):int(intergenic[i+1][1])]])
    finalwrite = []
    for i in range(len(extractintergenic)):
        if extractintergenic[i][1] == "":
            pass
        else:
            finalwrite.append(extractintergenic[i])
    with open(outfile, "w") as fastawrite:
        for i in range(len(finalwrite)):
            fastawrite.write(f">{finalwrite[i][0]}\n{finalwrite[i][1]}\n")
        fastawrite.close()
