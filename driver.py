"""
This is the driver program for Planted Motif Search by PSO
as implemeted in U. Srinivasulu Reddy, Michael Arock, A.V. Reddy, “A Particle Swarm Optimization Solution for Challenging Planted (l, d)-Motif Problem”, 2013 IEEE Symposium on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB)
"""
import pms_pso as p
def outFileName(string, index):
    return string[:index] + '_out' + string[index:]
def main():
	fileName=raw_input("Please Enter the File Name : ")
	seq=p.readInput(fileName)
	l=input("Enter the value of l ")
	d=input("Enter the value of d ")
	motif=list(p.pso(seq,l,d))
	ind=fileName.find('.')
	with open(outFileName(fileName,ind),'w') as outFile :
		for i in range(len(motif)):
       			outFile.write("Sequence "+str(i+1)+'\t'+str(motif[i])+'\n')

main()						
	
