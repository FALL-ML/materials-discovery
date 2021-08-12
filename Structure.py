# coding=utf8
'''
Created on 
@author: hebing
'''
import re
import numpy as np
import ase.spacegroup as spg
from scipy.spatial.ckdtree import cKDTree
#from pymatgen.core.periodic_table  import Element 
class Site(object):
    '''
    
    '''
    def __init__(self,sitelabel=[],Pos=None,elementsoccupy=None):
        self._Id=0
        self._SiteLabel=sitelabel
        self._SiteType=[]
        self._FracPostion=Pos
        self._Elements_occupy={}
        self._Elements_OxiValue={}
        self.SetElementsOccupy(elementsoccupy)
    def SetSiteType(self,sitetype):
        if sitetype not in self._SiteType:
            self._SiteType.append(sitetype)
    def GetSiteType(self):
        return self._SiteType
    def SetSiteLabel(self,sitelabel):
        if sitelabel != self._SiteLabel:
            self._SiteLabel=self._SiteLabel+sitelabel
    def SetPostion(self,pos=None):
        self._FracPostion=pos
    def SetElementsOccupy(self,eleocu):
        for key in eleocu.keys():
            if key in self._Elements_occupy:
                pass#self._Elements_occupy[key]=self._Elements_occupy[key]+eleocu[key]
            else:
                self._Elements_occupy[key]=eleocu[key]
    def SetElementsOxiValue(self,eleoxi):
        for key in eleoxi.keys():
                self._Elements_OxiValue[key]=eleoxi[key]
    def GetElementsOxiValue(self):
        return self._Elements_OxiValue
    def SetSiteId(self,siteid=0):
        self._Id=siteid
    def GetSiteId(self):
        return self._Id;
    def GetSiteLabel(self):
        return self._SiteLabel
    def GetPosition(self):
        return self._FracPostion
    def GetElementsOccupy(self ):
        return self._Elements_occupy
    def GetElements(self):
        return self._Elements_occupy.keys()
    def GetIronType(self):
        for (key,value) in self._Elements_OxiValue.items():
            if value<0:
                return -1
            elif value>0:
                return 1
            else:
                return 0
        return 0
class PeriodSite(Site):
    def __init__(self,parent=None,sitelabel=[],sitetype=None,Pos=None,elementsoccupy=None):
        super(PeriodSite,self).__init__(sitelabel,sitetype,Pos, elementsoccupy)
        self.__ParentStruct=parent
        self.__ExpandCell=np.array([[0,0,0],[0,0,0]])
    def SetParentStructure(self,parent):
        self.__ParentStruct=parent
    def GetExpandCell(self):
        return self.__ExpandCell
class Bond(object):
    def __init__(self,site0=None,site1=None,startpos=[0,0,0],endpos=[0,0,0]):   
        self._Site0=site0
        self._Site1=site1
        self._Positions = [startpos, endpos]

    def get_postions(self):
        return self._Positions


    def set_postions(self, value):
        self._Positions = value


    def del_postions(self):
        del self._Positions


    def get_site_0(self):
        return self._Site0

    def get_site_1(self):
        return self._Site1

    def set_site_0(self, value):
        self._Site0 = value

    def set_site_1(self, value):
        self._Site1 = value

    def del_site_0(self):
        del self._Site0

    def del_site_1(self):
        del self._Site1

    def SetSites(self,site0=None,site1=None):
        self._Site0=site0
        self._Site1=site1   
    Site0 = property(get_site_0, set_site_0, del_site_0, "Site0's docstring")
    Site1 = property(get_site_1, set_site_1, del_site_1, "Site1's docstring")
    Positions = property(get_postions, set_postions, del_postions, "Postions's docstring")
         
class Polyhedra(object):
    def __init__(self):
        self._VertexSites=[]
        self._CentreSite=None
        self._VertexVectors=[]

    def get_vertex_vectors(self):
        return self._VertexVectors


    def set_vertex_vectors(self, value):
        self._VertexVectors = value


    def del_vertex_vectors(self):
        del self._VertexVectors


    def get_centre_site(self):
        return self._CentreSite


    def set_centre_site(self, value):
        self._CentreSite = value


    def del_centre_site(self):
        del self._CentreSite


    def get_sites(self):
        return self._VertexSites


    def set_sites(self, value):
        self._VertexSites = value


    def del_sites(self):
        del self._VertexSites[:]

    def get_label(self):
        return self._CentreSite.GetSiteLabel()
    Label = property(get_label)
    VertexSites = property(get_sites, set_sites, del_sites, "Sites's docstring")
    CentreSite = property(get_centre_site, set_centre_site, del_centre_site, "CentreSites's docstring")
    VertexVectors = property(get_vertex_vectors, set_vertex_vectors, del_vertex_vectors, "VertexVectors's docstring")
class Structure(object):
    def __init__(self):
        self.__abc=[]
        self.__alphabetagama=[]
        self.__positions=[]
        self._initsites=[]
        self._Sites=[]
        self._SuperCellusites=[]
        self._SuperCellusitesFracPos=[]
        self._SuperCellusitesCartPos=[]
        self.__SuperCellusitesCellIndex=[]
        self.__spacegroupno=1
        self.__ABC=[]
        self.__equivviewsites=[] 
        self.__KDtree=None
 #       self.__ElementsSymbol={}
        self._Polyhedras={}
        self._Bonds={}
        self.__usitepos=[]
        self._OxidationTable={}
        self._atomsymbols=[]
    def get_atomsymbols(self):
        return self._atomsymbols
    def get_polyhedras(self):
        return self._Polyhedras
    def set_polyhedras(self, value):
        self._Polyhedras = value
    def del_polyhedras(self):
        del self._Polyhedras

    def get_bonds(self):
        return self._Bonds

    def set_bonds(self, value):
        self._Bonds = value

    def del_bonds(self):
        del self._Bonds

    def Setabc(self,abc):
        self.__abc=abc
    def Getabc(self):
        return self.__abc
    def SetABC(self,ABC=None):
        self.__ABC=ABC
    def GetABC(self):
        return self.__ABC
    def GetSuperCellusitesFracPos(self):
        return self._SuperCellusitesFracPos
    def GetSuperCellusites(self):
        return self._SuperCellusites
    def GetSuperCellusitesCartPos(self):
        return self._SuperCellusitesCartPos
    def ExpandCell(self):
        indexpos=[-2.0,-1.0,0.0,1.0,2.0]
        del self._SuperCellusites[:]
        del self._SuperCellusitesFracPos[:]
        del self._SuperCellusitesCartPos[:]
        for k in range(5):
            for j in range(5):
                for i in range(5):
                    for ste in self._Sites:
                        pos=ste.GetPosition()
                        stepos=[pos[0]+indexpos[i],pos[1]+indexpos[j],pos[2]+indexpos[k]]
                        self._SuperCellusitesFracPos.append(stepos.copy())
                        cartpos=self.FracPosToCartPos(stepos)
                        self._SuperCellusitesCartPos.append(cartpos.copy())
                        self._SuperCellusites.append(ste) 
        self.__KDtree=cKDTree(self._SuperCellusitesCartPos)
    def GetNeighbors(self,centre,r):  
        return self.__KDtree.query_ball_point(centre, r)
    def GetKNeighbors(self,centre,kn):  
        return self.__KDtree.query(centre, kn,p=2,n_jobs=-1)
    def FracPosToCartPos(self,FracPos):
        #pos=self.__ABC[0]*FracPos[0]+self.__ABC[1]*FracPos[1]+self.__ABC[2]*FracPos[2]
        pos=np.dot(FracPos,self.__ABC)
        return pos
    def CreateBonds(self,r):
        self._Bonds.clear()
        for ste in self._Sites:
            if ste.GetIronType()>0:
                    stepos=ste.GetPosition()
                    cartpos=self.FracPosToCartPos(stepos)
                   # neighbors=self.GetKNeighbors(cartpos, r=3)
                    (distance,neighbors)=self.GetKNeighbors(cartpos, kn=12)
                    if len(neighbors)>1:
                        polyhedra=Polyhedra()
                        polyhedra.CentreSite=ste
                        
                        for (i,index) in enumerate(neighbors[1:]):
                            site2=self._SuperCellusites[index]
                            fracpos2=self._SuperCellusitesFracPos[index]
                            if site2.GetIronType()<0 and distance[i+1]<r:
                                bond=Bond(ste,site2)
                                bond.Positions=[stepos,fracpos2]
                                if not ste.GetSiteId() in self._Bonds:
                                    self._Bonds[ste.GetSiteId()]=[]
                                self._Bonds[ste.GetSiteId()].append(bond)
                                polyhedra.VertexSites.append(site2)
                                polyhedra.VertexVectors.append(fracpos2-stepos)
                            else:
                                break
                        self._Polyhedras[ste.GetSiteId()]=polyhedra
    def CreateAllBonds(self):
        self.CreateBonds(r=3)
               
    def SetSpaceGroupno(self,spacegroupid=1):
        self.__spacegroupno=spacegroupid
    def GetSpaceGroupno(self):
        return self.__spacegroupno
    def SetAlphaBetaGama(self,abg):
        self.__alphabetagama=abg
    def GetAlphaBetaGama(self):    
        return self.__alphabetagama
    def AddSite(self,site=None):
        self._initsites.append(site)
    def GetSites(self):
        return self._initsites
    def GetUSites(self):
        return self._Sites
    def GetAseStructure(self,atoms=None):
        self.Setabc([atoms.info['_cell_length_a'],atoms.info['_cell_length_b'],atoms.info['_cell_length_c']])
        self.SetAlphaBetaGama([atoms.info['_cell_angle_alpha'],atoms.info['_cell_angle_beta'],atoms.info['_cell_angle_gamma']])
        if '_symmetry_int_tables_number' in atoms.info:
            self.SetSpaceGroupno( atoms.info['_symmetry_int_tables_number' ])
        elif '_space_group_it_number' in atoms.info:
            self.SetSpaceGroupno( atoms.info['_space_group_it_number' ])
        else:
            print('can not find group no')
            return
        index=0
        self.__ABC=atoms.cell
        postions=list(zip(atoms.info['_atom_site_fract_x'],atoms.info['_atom_site_fract_y'],atoms.info['_atom_site_fract_z']))
       
        
        no = None
        if '_space_group.it_number' in atoms.info:
            no = atoms.info['_space_group.it_number']
        elif '_space_group_it_number' in atoms.info:
            no = atoms.info['_space_group_it_number']
        elif '_symmetry_int_tables_number' in atoms.info:
            no = atoms.info['_symmetry_int_tables_number']
    
        symbolHM = None
        if '_space_group.Patterson_name_h-m' in atoms.info:
            symbolHM = atoms.info['_space_group.patterson_name_h-m']
        elif '_symmetry_space_group_name_h-m' in atoms.info:
            symbolHM = atoms.info['_symmetry_space_group_name_h-m']
        elif '_space_group_name_h-m_alt' in atoms.info:
            symbolHM = atoms.info['_space_group_name_h-m_alt']

        sitesym = None
        for name in ['_space_group_symop_operation_xyz',
                 '_space_group_symop.operation_xyz',
                 '_symmetry_equiv_pos_as_xyz']:
            if name in atoms.info:
                sitesym = atoms.info[name]
                break
            else:
                sitesym = None
 
        spgroup = 1
        if sitesym is not None:
            # subtrans = [(0.0, 0.0, 0.0)]
            # spgroup = spg.spacegroup.spacegroup_from_data(no=no, symbol=symbolHM, sitesym=sitesym, subtrans=subtrans)
            spgroup=spg.Spacegroup(self.GetSpaceGroupno())
        
#        spgroup=spg.Spacegroup(self.GetSpaceGroupno())
        kinds=[]
        ssites,kinds=spgroup.equivalent_sites(postions,onduplicates='keep')
        self.__positions=ssites
        if '_atom_type_oxidation_number' in atoms.info:
            for i,atomtype  in enumerate(atoms.info['_atom_type_symbol']):
                self._OxidationTable[atomtype]=int(atoms.info['_atom_type_oxidation_number'][i])
        atomsymbols=set([])
        for ssite in ssites:
            s=atoms.info['_atom_site_type_symbol'][kinds[index]]
            m = re.search(r'([A-Z][a-z]?)', s)
            symbol = m.group(0)
            atomsymbols.add(symbol)
            stelabel=atoms.info['_atom_site_label'][kinds[index]]
            eleoccupy={symbol:atoms.info['_atom_site_occupancy'][kinds[index]]}
            self.AddSite(Site(sitelabel=stelabel,Pos=ssite,elementsoccupy=eleoccupy))
            self._initsites[index].SetSiteType(s)
            self._initsites[index].SetSiteId(index)
            if len(self._OxidationTable) >0:
                eleoxi={symbol:self._OxidationTable[s]}
                self._initsites[index].SetElementsOxiValue(eleoxi)
            index=index+1
        self._atomsymbols=[i for i in atomsymbols]
        self.MergeDuplicateSite(0.001)
        self.ExpandCell()
        self.CreateAllBonds()
    def MergeDuplicateSite(self,symprec=0.001):
        dupsites=[False]*len(self.__positions)
        sitelabels=[]
        for index,sitepos in enumerate(self.__positions):
            if not dupsites[index]:
                self.__usitepos.append(self.__positions[index])
                #_site=Site(sitelabel=self._initsites[index].GetSiteLabel(),Pos=sitepos,elementsoccupy=self._initsites[index].GetElementsOccupy())
                self._Sites.append(self._initsites[index])
                sitelabels.append([self._initsites[index].GetSiteLabel()])
                self._Sites[-1].SetSiteId(len(self._Sites)-1)
                t = sitepos - self.__positions
                mask = np.all((abs(t) < symprec) |(abs(abs(t) - 1.0) < symprec), axis=1)
                if np.any(mask):
                    inds = np.argwhere(mask).tolist()
                    if len(inds)>1:
                        for dindex in inds:
                            dupsites[dindex[0]]=True
                            _occupy=self._initsites[dindex[0]].GetElementsOccupy()
                            _sitelabel=self._initsites[dindex[0]].GetSiteLabel()
                            _sitetype=self._initsites[dindex[0]].GetSiteType()
                            if index !=dindex[0] and (not _sitelabel in sitelabels[-1] ):
                                sitelabels[-1].append(_sitelabel)
                                for stype in _sitetype:
                                   # print(index)
                                    self._Sites[-1].SetSiteType(stype)
                                    m = re.search(r'([A-Z][a-z]?)', stype)
                                    symbol = m.group(0)
                                    if len(self._OxidationTable) >0:
                                        eleoxi={symbol:self._OxidationTable[stype]}
                                        self._Sites[-1].SetElementsOxiValue(eleoxi)
                                self._Sites[-1].SetElementsOccupy(_occupy)
                                self._Sites[-1].SetSiteLabel(_sitelabel)

    Bonds = property(get_bonds, set_bonds, del_bonds, "Bonds's docstring")
    Polyhedras = property(get_polyhedras, set_polyhedras, del_polyhedras, "Polyhedras's docstring")
