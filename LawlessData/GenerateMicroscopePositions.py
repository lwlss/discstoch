import math, numpy
from scipy import optimize


def makePosition(x,y,z,lab):
    string='''{
         "GRID_COL": 0,
         "DEVICES": [
            {
               "DEVICE": "ZStage",
               "AXES": 1,
               "Y": 0,
               "X": %s,
               "Z": 0
            },
            {
               "DEVICE": "XYStage",
               "AXES": 2,
               "Y": %s,
               "X": %s,
               "Z": 0
            }
         ],
         "PROPERTIES": {},
         "DEFAULT_Z_STAGE": "ZStage",
         "LABEL": "%s",
         "GRID_ROW": 0,
         "DEFAULT_XY_STAGE": "XYStage"
      }'''%(z,y,x,lab)
    return(string)

tlX,tlY,tlZ=-287411.87,203636.02,-15654.998
trX,trY,trZ=-368358.63,203533.81,-15791.636
brX,brY,brZ=-368755.7,159745.53,-15914.048
blX,blY,blZ=-287997.42,158808.89,-15708.12

nC,nR=10,6

def zval(x,y,a,b,d):
    return(a*x+b*y+d)

def zobj(xvec):
    [a,b,d]=xvec
    tl=zval(tlX,tlY,a,b,d)-tlZ
    tr=zval(trX,trY,a,b,d)-trZ
    br=zval(brX,brY,a,b,d)-brZ
    bl=zval(blX,blY,a,b,d)-blZ
    return(math.sqrt(tl**2+tr**2+br**2+bl**2))

meanz=(tlZ+trZ+brZ+blZ)/4.0
test=optimize.fmin_l_bfgs_b(zobj,numpy.array([0,0,meanz]),approx_grad=True,bounds=[(-10,10),(-10,10),(0.5*meanz,2.0*meanz)],maxfun=1000,epsilon=0.5)
    

#zave=(float(tlz)+float(brz))/2.0
#[trX,trY,blX,blY]=getCorners([tlX,tlY,brX,brY],12,8)

outfile=open("PythonMicroPositions.pos","w")
outfile.write('''{
   "VERSION": 3,
   "ID": "Micro-Manager XY-position list",
   "POSITIONS": [
   ''')
count=0
for R in xrange(0,nR):
    if R%2==0:
        cfrom,cto,delta=0,nC,1
    else:
        cfrom,cto,delta=nC-1,-1,-1
    for C in xrange(cfrom,cto,delta):
        lab="R%02dC%02d"%(R+1,C+1)
        cfrac=(float(C)/float(nC-1))
        rfrac=(float(R)/float(nR-1))
        x=tlX*(1-cfrac)*(1-rfrac)+trX*cfrac*(1-rfrac)+blX*(1-cfrac)*rfrac+brX*cfrac*rfrac
        y=tlY*(1-cfrac)*(1-rfrac)+trY*cfrac*(1-rfrac)+blY*(1-cfrac)*rfrac+brY*cfrac*rfrac
        z=tlZ*(1-cfrac)*(1-rfrac)+trZ*cfrac*(1-rfrac)+blZ*(1-cfrac)*rfrac+brZ*cfrac*rfrac
        outfile.write(makePosition(x,y,z,lab))
        print(lab,str(x),str(y),str(z))
        if(count<((nC)*(nR))):
            outfile.write(",")
            count+=1
outfile.write('''
]
}''')
outfile.close()



