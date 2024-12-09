path = 'aligned-genes/'
files = ['BDNF', 'FOXP2', 'MBP', 'OPN1SW', 'TBXT']

for file in files:
    sequence_names = []
    sequences = {}
    with open (path + file + '.fasta') as f:
        line = f.readline().strip()
        while line != '':
            if '>' in line:
                line = line.split(' ')
                if 'PREDICTED:' in line: sequence_names += ['>' + ' '.join(line[3:5])]
                else: sequence_names += ['>' + ' '.join(line[2:4])]
                sequences[sequence_names[-1]] = ''
            else:
                sequences[sequence_names[-1]] += line
            line = f.readline().strip()
        
    with open(path + file + '.fasta', 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')
            f.write(sequences[seq] + '\n')
