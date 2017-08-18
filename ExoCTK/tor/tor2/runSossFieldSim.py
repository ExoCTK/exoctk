"""
This file is meant to be used as a module for the ExoCTK website i.e. 
app_exoctk.py. It is part of the ExoCTK package. It is designed to test
the sossFieldSim.py file which runs the simulation. It is meant to run 
from the terminal. Usage of this file is described in sossFieldSim.py 
documentation.

Authors:
	Rafia Bushra, University of Arizona
	
	David Lafreniere, University de Montreal

"""

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
