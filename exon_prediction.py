import numpy as np

path = 'genes/'
alignedpath = 'realigned-genes/'
testpath = 'test-realigned-genes/'
files = ['BDNF', 'FOXP2', 'MBP', 'OPN1SW', 'TBXT', 'RBFOX1']

for file in files:
    sequence_names = []
    sequences = {}
    with open (path + file + '_refseq_transcript' + '.fasta') as f:
        line = f.readline().strip()
        while True:
            if not line: break
            line = line.strip()
            if '>' in line:
                line = line.split(' ')
                if 'PREDICTED:' in line: sequence_names += ['>' + '-'.join(line[2:4])]
                else: sequence_names += ['>' + '-'.join(line[1:3])]
                sequences[sequence_names[-1]] = ''
            else:
                sequences[sequence_names[-1]] += line
            line = f.readline()

    minLen = 100000
    minSeqName = ""
    for seq in sequences:
        if len(sequences[seq]) < minLen:
            minSeqName = seq
            minLen = len(sequences[seq])

    aligned_sequences = {}
    with open(alignedpath + file + '.fasta') as f:
        line = f.readline().strip()
        name = ""
        while line != '':
            if '>' in line: 
                name = line
            else:
                aligned_sequences[name] = line
            line = f.readline().strip()
    
    with open(testpath + file + '.fasta', 'w') as f:
        minSeq = aligned_sequences[minSeqName]
        k = 2
        for seq in sequences:
            currSeq = ""
            f.write(seq + '\n')
            for pos in range(len(minSeq)):
                if minSeq[pos : pos + k] != ('-' * k): currSeq += aligned_sequences[seq][pos]
            f.write(currSeq + "\n")

def find_genetic_distances(infile, outfile = 'genetic-distances.txt'):
    sequenceNames = []
    sequences = {}

    with open(infile) as f: # Read in input fna file
        line = f.readline()
        while line != '':
            if line[0] == '>': sequenceNames.append(line[1:].strip()) # Save sequence names in list
            else: sequences[sequenceNames[-1]] = line.strip() # Save sequences in dictionary
            line = f.readline()

    distances = np.zeros((len(sequenceNames), len(sequenceNames))) # Initialize distances array
    distanceStrings = np.empty_like(distances, dtype = object) # Copy of distances array as strings

    numSequences = len(sequenceNames)

    for i in range(numSequences):
        for j in range(numSequences):
            if i == j:
                distances[i, j] = 0.0 # Distance from a sequence to itself is 0
            else:
                seq1 = sequences[sequenceNames[i]]
                seq2 = sequences[sequenceNames[j]]
                seqLen = len(seq1)

                matches = 0 # Count matches and save distance
                for k in range(seqLen):
                    if seq1[k] == seq2[k]: matches += 1
                
                distances[i, j] = 1 - (matches / seqLen) 
            distanceStrings[i, j] = str(distances[i, j]) # Write to distancesStrings

    with open(outfile, 'w') as f:
        header = '\t'
        header += '\t'.join(sequenceNames) + '\n'
        f.write(header)

        for i in range(numSequences): # Write out sequence name and genetic distance info to file
            lineToWrite = sequenceNames[i] + '\t'
            lineToWrite += '\t'.join(distanceStrings[i, :]) + '\n'
            f.write(lineToWrite)
    return distances

find_genetic_distances('test-realigned-genes/BDNF.fasta')