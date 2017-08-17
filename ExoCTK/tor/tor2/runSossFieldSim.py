from sys import argv
from sossFieldSim import *


if len(argv) == 3:
	sossFieldSim(argv[1],argv[2])

elif len(argv) == 4:
	sossFieldSim(argv[1],argv[2],argv[3])

elif len(argv) == 5:
	 bc  = float(argv[4].split(','))
	 sossFieldSim(argv[1], argv[2], argv[3], bincomp=bc)

else:
	print("Needs 2, 3 or 4 arguments")
