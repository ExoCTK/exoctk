# !path = '+/home/david/idl/mylib/:' + !path
# !path = expand_path(!path)

from sys import argv

from . import sossFieldSim

if len(argv) == 3:
	sossFieldSim.sossFieldSim(argv[1],argv[2])
elif len(argv) == 4:
   
	 bc  = float(argv[3].split(','))
   
	 sossFieldSim.sossFieldSim(argv[1], argv[2], bincomp=bc)
else:
	print("Needs 2 or 3 arguments")
