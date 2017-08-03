# This takes a file and changes order of the columns
# example
# python process.py Desktop/LN_CD8+CD69+.txt Desktop/processed.txt
import sys
file_1 = sys.argv[1]
file_1_new = sys.argv[2]

with open(file_1, 'r') as startingfile:
	with open(file_1_new), 'w') as endingfile:
		endingfile.write("vMaxResolved.jMaxResolved\tcount\tvMaxResolved\tjMaxResolved\taminoAcid\tnucleotide\n")
		startingfile.readline() #read the first line of the file and don't do anything with it
		for row in startingfile:
			data = row.split()
			output.write(data[4] + "\t" + data[0] + "\t" + data[4] + "\t" + data[6] + "\t" + data[3] + "\t" + data[2] + "\n")
