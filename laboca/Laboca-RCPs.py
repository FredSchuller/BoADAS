####################################
def flagRCP(rcpFile):
    
    try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile))
    except IOError:
            self.__MessHand.error("could not open file %s"%(rcpFile))
            return

    param = f.readlines()
    f.close()

    channels = []
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
                tmp = string.split(param[i])
                num = string.atoi(tmp[0])
                channels.append(num)

    missing = []
    
    for i in range(data.BolometerArray.NChannels):
        b = i+1
        found = 0
        for j in range(len(channels)):
            if channels[j] == b:
                found = 1

        
        if found == 0:
            missing.append(b)

    flagC(missing)

def getLabocaRCP(mjdref):
	
	try:
 		f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-rcps.dat'))
	except IOError:
		self.__MessHand.error("could not open Laboca-rcps.dat")	
		return

	param = f.readlines()	
	f.close()

	mjdstart = []
	mjdstop = []
	rcp=[]
	for i in range(len(param)):
 		if param[i][0] not in ['!','#','n']:    # skip comments
                	tmp = string.split(param[i])
                	mjdstart.append(string.atof(tmp[0]))
                	mjdstop.append(string.atof(tmp[1]))
			sr=tmp[2]
                	rcp.append(sr)

	mjdstart = array(mjdstart)
    	mjdstop =  array(mjdstop)
    	
	rcpname= -1

	for i in range(len(mjdstart)):
		low = mjdstart[i]
		high = mjdstop[i]
		if low <= mjdref and mjdref < high:
			rcpname = rcp[i]

	return rcpname


def getLabocaCross(mjdref):
	
	try:
 		f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-cross.dat'))
	except IOError:
		self.__MessHand.error("could not open Laboca-cross.dat")	
		return

	param = f.readlines()	
	f.close()

	mjdstart = []
	mjdstop = []
	cr=[]
	for i in range(len(param)):
 		if param[i][0] not in ['!','#','n']:    # skip comments
                	tmp = string.split(param[i])
                	mjdstart.append(string.atof(tmp[0]))
                	mjdstop.append(string.atof(tmp[1]))
			sr=tmp[2]
                	cr.append(sr)

	mjdstart = array(mjdstart)
    	mjdstop =  array(mjdstop)
    	
	cross= []
        crossreturn = []
	for i in range(len(mjdstart)):
		low = mjdstart[i]
		high = mjdstop[i]
		if low <= mjdref and mjdref < high:
			cross = cr[i]


        tmp=cross.split(',')
        for a in tmp:
            crossreturn.append(int(a))

	return crossreturn


