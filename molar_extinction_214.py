#!/usr/bin/python
'''Calculate molar extinction coefficient of a peptide
at 214 nm based on primary AA sequence, created by David Langelaan'''

#Insert the sequence of your peptide here...no spaces
sequence = 'GVFVHNSAGCISGDSLITLA'


#Don't really need an output file.
#Outputfilename = 'Spitz_214nm_extinction.txt'

#individual molar extinction coefficients
#Proline N-terminal is 30 however, will need an if statement for that?

Residues = {'P': 2675,'H': 5125,'F': 5200,'Y': 5375,'W': 29050,'M': 980,'R': 102,'N': 136,'Q': 142,
	    'C': 225,'G': 21,'A': 32,'S': 34,'K': 41,'T': 41,'V': 43,'I': 45,'L': 45,'D': 58,'E': 78}
	    
#peptide bond molar extinction coefficient
peptide_bond = 923

#Let e = extinction coefficient

# EQUATION:
#e  = (923)*(n-1)+(frequency of AA i)*(e value of AA i)

# where (n-1) is the number of peptide bonds (i.e. number of letters in the AA sequence string minus 1)
# will need a check for N-terminal proline

absorbances = []
for index, residue in enumerate(sequence):
    if index == 0 and residue == 'P':
        absorbances.append(30)
    else:
        absorbances.append(Residues[residue])

print peptide_bond*(len(sequence)-1) + sum(absorbances)
    
	    

