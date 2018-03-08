import os

def PlanetFlux(source,time,date,beam,freq):
    sabocadir = os.getenv('BOA_HOME_SABOCA')
    script= sabocadir+'/saboca-planet-flux.scr'
    input = [script,str(source),str(time),str(date),str(beam),str(freq)]

    command = " ".join(input)
    os.system(command)

    f = file(os.path.join(os.getenv('HOME'),'saboca-planet-flux.dat'))
    rl = f.readlines()
    f.close()

    result=[]
    for l in rl:
        tmp = l.split()
        result.append(float(tmp[0]))

        flux = result[0]

    return flux


def getAstroDate(data):
    months = ['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']
    dateObs = data.ScanParam.DateObs
    month = int(dateObs[5:7])
    year = dateObs[:4]
    day = dateObs[8:10]
    dateAstro = day + '-' + months[month-1] + '-' + year
    timeAstro = dateObs[11:] + '.0'
    return timeAstro,dateAstro
    
