#!/usr/bin/env python
import sys
import string
#takes a list of pdbs with characters of the chains you want to keep
# example 1uad.pdb A,C

if __name__ == '__main__':
    pdbchainlist = open(sys.argv[1])

    for eachline in pdbchainlist:
        if len(eachline) < 2:
            print "Line is empty.  Moving on..."
            continue
        splitline = eachline.split()
        # make the file name
        filename = splitline[0]
        output_filename= filename[:-4] + "_chain_"
        inputpdb = open( filename , "r" )
        pdb_str = inputpdb.read()
        commachains = splitline[1]
        eachchain=commachains.split(',')
        keptlines = []
        #itterate through chains you want to keep
        for chain in eachchain:
        #itterate through the pdb
            for l in pdb_str.split("\n"):
                if l[21:22] == chain:
                    #keep any MSE lines
                    if l[:6] == "HETATM" and l[17:20] == "MSE":
                        #fixes missing occupancy
                        line = l[:56]+"1.00"+l[60:]
                        keptlines.append(line)
                    #other residues
                    elif l[:4] == "ATOM" or l[:3] == "TER":
                        #fixes missing occupancy
                        line = l[:56]+"1.00"+l[60:]
                        keptlines.append(line)
        #finish itterating and building 

        #name output file
        for substr in eachchain:
            output_filename = output_filename + substr
     
        output_filename = output_filename + ".pdb"

        print "Writing file: ",output_filename

        #now open and write to said file
        newpdb = open( output_filename , "w" )
      
        #write the mutation(s) to the file
        for line in keptlines:
            newpdb.write(line + "\n")
        #close the file when done
        newpdb.close()
        #completed each line
        continue
    print "Completed!"



        
