#python script to creat a list of cases, where by the cases are the females
#


import sys

fam = open(sys.argv[1], 'r')
pheno = open(sys.argv[2], 'w')

fam_a = fam.readlines()

for i in fam_a:
	line=i.split()
	if line[4] == '2':#check fam file 5th colums is sex usually
		pheno.write('{}\n'.format(line[0]))
