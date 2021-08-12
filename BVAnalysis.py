'''
Created on 2017

@author: hebing
'''
import os
import numpy as np
from scipy.special import erfc 
from scipy.ndimage import measurements
import scipy.ndimage
import math
from struct import pack
def loadbvparam(filename):
    with open(filename,'r') as fp:
        lines=fp.readlines()
    bvmpara={}
    for line in lines:
        varstr=line.split('\t')
        key=varstr[0]+varstr[1]+varstr[2]+varstr[3]
        if key in bvmpara:
            bvmpara[key].append([float(varstr[4]),float(varstr[5])])
        else:
            bvmpara[key]=[[float(varstr[4]),float(varstr[5])]]
    return bvmpara
  
def loadBVSEparam(filename):
    with open(filename,'r') as fp:
        lines=fp.readlines()
    BVSEparam={}
    for line in lines[1:]:
        line=line.strip('\n')
        varstr=line.split('\t')
        key=varstr[0]+varstr[1]+varstr[2]+varstr[3]
        if key in BVSEparam:
            BVSEparam[key].append([float(varstr[4]),float(varstr[5]),float(varstr[6]),float(varstr[7]),float(varstr[8]),float(varstr[9])])
        else:
            BVSEparam[key]=[[float(varstr[4]),float(varstr[5]),float(varstr[6]),float(varstr[7]),float(varstr[8]),float(varstr[9])]]
    return BVSEparam

def loadElementsparam(filename):
    with open(filename,'r') as fp:
        lines=fp.readlines()
    Elementsparam={}
    for line in lines:
        line=line.strip('\n')
        varstr=line.split('\t')
        key=varstr[1]+varstr[2]
        Elementsparam[key]=[float(varstr[0]),float(varstr[3]),float(varstr[4]),float(varstr[5]),float(varstr[6]),float(varstr[7]),float(varstr[8]),float(varstr[9])]
    return Elementsparam    
class BVAnalysis(object):
    '''
    for bond valence method analyze ion channel
    '''
    def __init__(self, struc=None,size=[50,50,50]):
        '''
        Constructor
        '''
        self._Struct=struc
        self._Size=[50,50,50]
        self.__Data={}
        self.__Max={}
        self.__Min={}
        self.Ea={}
        self._MoveIon='Li'
        self._poss=None
        self.ValenceOfMoveIon=1
        self.stop=False
        module_dir = os.path.dirname(os.path.abspath(__file__))

        self._BVparam= loadbvparam(os.path.join(module_dir, 'bvmparam.dat'))
        self._BVSEparam= loadBVSEparam(os.path.join(module_dir, 'bvse.dat'))
        self._Elementsparam=loadElementsparam(os.path.join(module_dir, 'elements.dat'))

    def get_data(self):
        return self.__Data


    def get_max(self):
        return self.__Max


    def get_min(self):
        return self.__Min


    def set_data(self, value):
        self.__Data = value


    def set_max(self, value):
        self.__Max = value


    def set_min(self, value):
        self.__Min = value


    def del_data(self):
        del self.__Data


    def del_max(self):
        del self.__Max


    def del_min(self):
        del self.__Min


    @property
    def BVSMin(self):
        return self.__Min['BVS']
    @property
    def BVSEMin(self):
        return self.__Min['BVSE']
    @property
    def BVELMin(self):
        return self.__Min['BVEL']
    @property
    def BVSMax(self):
        return self.__Max['BVS']
    @property
    def BVSEMax(self):
        return self.__Max['BVSE']
    @property
    def BVELMax(self):
        return self.__Max['BVEL']
    @property
    def BVSEEa(self):
        return self.Ea['BVSE']
    @property
    def BVELEa(self):
        return self.Ea['BVEL']
    def SetStructure(self,struc):
        self._Struct=struc
    def SetMoveIon(self,moveion):
        self._MoveIon=moveion
   
    def SetResolution(self,size):
        self._Size=size
    def SetLengthResolution(self,reslen=0.1):
        abc=self._Struct.Getabc()
        self._reslen=reslen
        self._Size[0]=int(abc[0]/reslen)+1
        self._Size[1]=int(abc[1]/reslen)+1
        self._Size[2]=int(abc[2]/reslen)+1
    def SaveBVSEData(self,filename):
        self.SaveBinaryData(filename+'_BVSE.pgrid', self.__Data['BVSE'],'BVSE')
        np.save(filename+"_BVSE.npy",self.__Data['BVSE'])
    def SaveBVELData(self,filename):
        self.SaveBinaryData(filename+'_BVEL.pgrid', self.__Data['BVEL'],'BVEL')
    def SaveBVSData(self,filename):
        self.SaveBinaryData(filename+'_BVS.pgrid', self.__Data['BVS'],'BVS')
    def SaveData(self,filename,data,datatype=' '):
        fp=open(filename,'w')
        fp.write(datatype+'\n')
        alpha=self._Struct.GetAlphaBetaGama()
        abc=self._Struct.Getabc()
        fp.write(str(abc[0])+' '+str(abc[1])+' '+str(abc[2])+' ')
        fp.write(str(alpha[0])+' '+str(alpha[1])+' '+str(alpha[2])+'\n')
        fp.write(str(self._Size[0])+' '+str(self._Size[1])+' '+str(self._Size[2])+'\n ')
        for i in range(self._Size[0]):
            for j in range(self._Size[1]):
                for k in range(self._Size[2]):
                    fp.write(str(data[i][j][k])+'\n')
        fp.close()
    def SaveBinaryData(self,filename,data,datatype=' '):
        latticeparam=self._Struct.Getabc()+self._Struct.GetAlphaBetaGama()
        with open(filename, 'wb') as outf:
            outf.write(pack('4i', 3, 0, 0, 0))  # file format version
            outf.write(pack('80s', datatype.ljust(80).encode('utf-8')))  # Title 80 characters
            outf.write(pack('i', 1))  # gType 0 for *.ggrid, 1 for *.pgrid
            outf.write(pack('i', 0))  # fType 0
            outf.write(pack('i', 1))  # nVal dummy
            outf.write(pack('i', 3))  # dimension
            outf.write(pack('3i', *self._Size))  # numbers of voxels
            outf.write(pack('i', data.size))  # Total number of voxels
            outf.write(pack('6f', *latticeparam))  # lattice parameters
            for k in range(self._Size[2]):
                for j in range(self._Size[1]):
                    for i in range(self._Size[0]):
                        outf.write(pack('f', data[i][j][k]))
    def ReadData(self,filename,data,datatype=' '):
        fp=open(filename,'r')
        datatype=fp.readline()
        alpha=self._Struct.GetAlphaBetaGama()
        abc=self._Struct.Getabc()
        abc=fp.readline().split()
        alpha=fp.readline().split()
        strs=fp.readline().split()
        size=[ int(i) for i in strs]
        if datatype=='bvs':
            self.__Data['BVS']=np.zeros(size)
            data=self.__Data['BVS']
        else:
            self.__Data['BVSE']=np.zeros(size)
            data=self.__Data['BVSE']
        for i in range(size[2]):
            for j in range(size[1]):
                for k in range(size[0]):
                    strs=fp.readline().split()
                    data[k][j][i]=int(strs[0])
                 
        fp.close()
    def GetBVSData(self):
        return self.__Data['BVS']
    def GetBVSEData(self):
        return self.__Data['BVSE']

    def GetPosSet(self):
        return self._poss
    def CaluGlobalIndex(self,):
        self.stop=False
        self._Size.reverse()
        centrepos=np.zeros(self._Size+[3])
        self._Data=np.zeros(self._Size,dtype=np.double)
        for k in range(self._Size[0]):
            for j in range(self._Size[1]):
                for i in range(self._Size[2]):
                    centrepos[k][j][i][2]=k/(self._Size[0]-1.0)
                    centrepos[k][j][i][1]=j/(self._Size[1]-1.0)
                    centrepos[k][j][i][0]=i/(self._Size[2]-1.0)
 
        self._poss=self._Struct.FracPosToCartPos(centrepos)
        (distance,neighborsindex)=self._Struct.GetKNeighbors(self._poss, kn=8)
              
        site2atoms=None
        for k in range(self._Size[0]):
 
            for j in range(self._Size[1]):
                for i in range(self._Size[2]):
                    if self.stop:
                        return 
                    for dindex,index in enumerate(neighborsindex[k][j][i]):
                        site2=self._Struct.GetSuperCellusites()[index]
                        if site2.GetIronType()*self.ValenceOfMoveIon<0:
                            for key in site2.GetElementsOccupy():
                                if key!='Vac':
                                    site2atoms=key
                            if self.ValenceOfMoveIon>0:
                                key="".join([self._MoveIon,str(self.ValenceOfMoveIon),site2atoms,str(site2.GetElementsOxiValue()[site2atoms])])
                            else:
                                key="".join([site2atoms,str(site2.GetElementsOxiValue()[site2atoms]),self._MoveIon,str(self.ValenceOfMoveIon)])
                            if key in self._BVparam:
                                (r,b)=self._BVparam[key][0]
                                bv=np.exp((r-distance[k][j][i][dindex])/b)
                                self.__Data['BVS'][k][j][i]=self.__Data['BVS'][k][j][i]+bv
                            else:
                                print('Not Find bvs param'+key)
                        else:
                            pass
                            #print(site2.GetSiteLabel())                            
    def CaluBVSE(self,ProgressCall):
        if len(self._Struct._OxidationTable)==0:
            raise Exception("can't caculation bvs and Ea without atom oxi info !")
        self.stop=False
        centrepos=np.zeros(self._Size+[3])
        self.__Data['BVS']=np.zeros(self._Size,dtype=np.double)
        self.__Data['BVSE']=np.zeros(self._Size,dtype=np.double)
        self.__Data['BVEL']=np.zeros(self._Size,dtype=np.double)
        for k in range(self._Size[2]):
            for j in range(self._Size[1]):
                for i in range(self._Size[0]):
                    centrepos[i][j][k][2]=k/(self._Size[2]-1.0)
                    centrepos[i][j][k][1]=j/(self._Size[1]-1.0)
                    centrepos[i][j][k][0]=i/(self._Size[0]-1.0)
        if ProgressCall:
            ProgressCall(10)
        self._poss=self._Struct.FracPosToCartPos(centrepos)
                        
        atomsq={}
        for atom in self._Struct._atomsymbols:
            atomsq[atom]=0
            for site in self._Struct._Sites:
                if atom in site.GetElements():
                    key=atom+str(site.GetElementsOxiValue()[atom])
                    atomsq[atom]=atomsq[atom]+site.GetElementsOxiValue()[atom]*site.GetElementsOccupy()[atom]/math.sqrt(self._Elementsparam[key][3])
        qsumanion=0
        qsumcation=0
        for (atom,value) in atomsq.items():
            if value >0:
                qsumcation=qsumcation+value
            elif value<0:
                qsumanion=qsumanion+value
            else:
                qsumanion=0
                qsumcation=0
        qsum=0.0
        if self.ValenceOfMoveIon>0:
            qsum=-qsumanion/qsumcation 
        else:
            qsum=-qsumcation/qsumanion          
        key1=self._MoveIon+str(self.ValenceOfMoveIon)  
        qm1=self.ValenceOfMoveIon/math.sqrt(self._Elementsparam[key1][3])
        rm1=self._Elementsparam[key1][6]
        for i in range(self._Size[0]):
            (distance,neighborsindex)=self._Struct.GetKNeighbors(self._poss[i], kn=128)
            if ProgressCall:
                ProgressCall(10+i*90/(self._Size[0]-1))
            for j in range(self._Size[1]):
                for k in range(self._Size[2]):
                    if self.stop:
                        return False
                    bvsdata=0.0
                    data=0.0
                    cdata=0.0
                    bvelcdata=0.0
                    N=0
                    D0=0.0
                    Rcutoff=10.0
                    alpha=0.0
                    occupyvalue=0.0
                    for dindex,index in enumerate(neighborsindex[j][k]):
                        site2=self._Struct._SuperCellusites[index]
                        if site2.GetIronType()*self.ValenceOfMoveIon<0:
                            for (ssymbol,occupyvalue) in site2.GetElementsOccupy().items():
                                if ssymbol!='Vac':
                                    site2oxi=site2.GetElementsOxiValue()[ssymbol]
                                    if self.ValenceOfMoveIon>0:
                                        key="".join([self._MoveIon,str(self.ValenceOfMoveIon),ssymbol,str(site2oxi)])
                                    else:
                                        key="".join([ssymbol,str(site2oxi),self._MoveIon,str(self.ValenceOfMoveIon)])
                                    if (key in self._BVSEparam) and (key in self._BVparam ):
                                        (Nc,r0,Rcutoff,D0,Rmin,alpha)=self._BVSEparam[key][0]
                                        (r,b)=self._BVparam[key][0]
                                        Rcutoff=10.0
                                        if distance[j][k][dindex]<=Rcutoff:
                                            #bv=np.exp((r0-distance[k][j][i][dindex])*alpha)
                                            bvs=np.exp((r-distance[j][k][dindex])/b)
                                            smin=np.exp((Rmin-distance[j][k][dindex])*alpha)
                                            en_s=np.exp((Rmin-Rcutoff)*alpha)
                                            bvsdata=bvsdata+bvs #*occupyvalue-((en_s-1)**2-1)
                                            data=data+((smin-1)**2-1)*occupyvalue
                                    else:
                                        if (key not in self._BVSEparam) :
                                            raise Exception("bvse {0}  param can't find!!!!".format(key))
                                            ProgressCall(0)
                                        else:
                                            raise Exception("bvs {0} param can't find!!!!".format(key))
                                            ProgressCall(0)
                        else:
                            for  (ssymbol,occupyvalue) in site2.GetElementsOccupy().items():
                                if ssymbol!='Vac':
                                    site2oxi=site2.GetElementsOxiValue()[ssymbol]
                                    if ssymbol !=self._MoveIon:
                                        key=ssymbol+str(site2oxi)   
                                        rm2=self._Elementsparam[key][6] 
                                        qm2=site2oxi/math.sqrt(self._Elementsparam[key][3])
                                        rm1m2=distance[j][k][dindex]
                                        f=0.74
                                        if rm1m2>rm2:
                                            if rm1m2<Rcutoff:
                                                cdata=cdata+occupyvalue*qm1*qm2/rm1m2*erfc(rm1m2/(f*(rm1+rm2)))*qsum
                                                bvelcdata=bvelcdata+occupyvalue*qm1*qm2/rm1m2*(erfc(rm1m2/(f*(rm1+rm2)))-erfc(Rcutoff/(f*(rm1+rm2))))
                                        else:
                                            cdata=300
                                            bvelcdata=300
                    self.__Data['BVSE'][i][j][k]=(0.5*D0*(data)+14.4*cdata)
                    self.__Data['BVS'][i][j][k]=bvsdata
                    self.__Data['BVEL'][i][j][k]=(D0*(data)+14.4*bvelcdata)
        self.__Data['BVSE']=self.__Data['BVSE']#/self.ValenceOfMoveIon
        self.__Data['BVEL']=self.__Data['BVEL']#/self.ValenceOfMoveIon
        for key,data in self.__Data.items():
            self.__Max[key]=np.max(data)
            self.__Min[key]=np.min(data)
            ashape=self.__Data['BVSE'].shape
            self.ndata=np.zeros((2*ashape[0],2*ashape[1],2*ashape[2]))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        self.ndata[i*ashape[0]:(i+1)*ashape[0], \
                                                 j*ashape[1]:(j+1)*ashape[1], \
                                                 k*ashape[2]:(k+1)*ashape[2]]=\
                                                data-self.__Min[key]
            self.Ea[key]=self.GetEa(self.ndata)

        #print('BVSE Ea is {0},\nBVEL Ea is {1}'.format(self.BVSEEa,self.BVELEa))
        return True
    def  GetEa(self,data):
        return [self.GetConnectEa(data,1),\
                self.GetConnectEa(data,2),\
                self.GetConnectEa(data,3)]

    def  GetConnectEa(self,data,connectnumber):             
            minvalue=0.0
            maxvalue=20.0
            while  (maxvalue-minvalue)>0.01:
                testvalue=(maxvalue+minvalue)/2.0
                testdata=data<testvalue
                number=self.CaluRegionConnectivity(testdata)
                if number>=connectnumber:
                    maxvalue=testvalue
                else:
                    minvalue=testvalue
            if maxvalue >=20.0:
                return 20000
            else:
                return maxvalue

    def CaluRegionConnectivity(self,region):
        struct = scipy.ndimage.generate_binary_structure(3,3)
        lw,num=measurements.label(region,structure=struct)
        sliced = measurements.find_objects(lw)
        connectnumber=np.zeros(num,dtype=np.int)
        if(num > 0):
            for index,slic in enumerate(sliced):
                for i in range(len(slic)):
                    if ((slic[i].start==0) and slic[i].stop==region.shape[i]):
                        connectnumber[index]=connectnumber[index]+1
        return np.max(connectnumber)                
        
    Data = property(get_data, set_data, del_data, "Data's docstring")
    Max = property(get_max, set_max, del_max, "Max's docstring")
    Min = property(get_min, set_min, del_min, "Min's docstring")
