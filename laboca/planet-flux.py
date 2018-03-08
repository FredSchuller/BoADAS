import os

def PlanetFlux(source,time,date,beam,freq):

    script= os.getenv('BOA_HOME_LABOCA')+'/planet-flux.scr'
    input = [script,str(source),str(time),str(date),str(beam),str(freq)]

    command = " ".join(input)
    os.system(command)


    f = file(os.path.join(os.getenv('HOME'),'planet-flux.dat'))
    rl = f.readlines()
    f.close()

    result=[]
    for l in rl:
        tmp = l.split()
        result.append(float(tmp[0]))

        flux = result[0]

    return flux


def getAstroDate(data):
    monthes = ['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']
    dateObs = data.ScanParam.DateObs
    month = int(dateObs[5:7])
    year = dateObs[:4]
    day = dateObs[8:10]
    dateAstro = day+'-'+monthes[month-1]+'-'+year
    timeAstro = dateObs[11:]+'.0'
    return timeAstro,dateAstro
    
