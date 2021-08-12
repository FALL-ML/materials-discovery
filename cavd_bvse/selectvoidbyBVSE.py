'''
__author__ = "Penghui Mi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Penghui Mi"
__email__ = "phmim@shu.edu.cn"
__date__ = "Mar 19, 2019"
'''
# 通过BVSE能量筛选掉能量值较高的间隙点
import numpy as np
import math
from pymatgen.core.sites import PeriodicSite
from scipy.interpolate import RegularGridInterpolator
from cavd_bvse.mergecluster import Void, Channel


class SelectVoidByBVSE(object):
    def __init__(self,struc,energy):
        self._energy = energy
        self._struc = struc
        self._threshold = 2000
        self._voids = []
        self._channels = []
        self._selectedvoids = []
        self._selectedchannels = []
        self.voidid_label = {}
        self.nonequalchannels = {}

    def setpara(self, voidslist,channelslist,bvse_Ea):
        # 读取BVSE方法计算出的能量场
        self._voids = voidslist
        self._channels = channelslist
        self._threshold = max([value for value in bvse_Ea])+np.amin(self._energy)

    @property
    def voids(self):
        # 返回结果为一个字典
        return self._selectedvoids

    @property
    def channels(self):
        # 返回结果为一个字典
        return self._selectedchannels

    def get_dis(self, p1, p2):
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis


    def calpointenergy(self, point):
        # 计算能量场中一点的能量
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z), self._energy)
        point = [math.modf(i)[0] for i in point]
        for i in range(len(point)):
            if point[i] < 0.0:
                point[i] += 1
        return interpolating_function(np.array(point))[0]

    def selectvoid(self):
        for void in self._voids:
            void.energy = self.calpointenergy(void.coord)
            if void.energy < self._threshold:
                self._selectedvoids.append(void)
                self.voidid_label[void.id] = void
        for channel in self._channels:
            if channel.start in self.voidid_label.keys() and channel.end in self.voidid_label.keys():
                coord1 = self.voidid_label[channel.start].coord
                coord2 = self.voidid_label[channel.end].coord
                channel.dist = self.get_dis(coord1, coord2)

                channel_energy = self.calpointenergy(channel.coord)
                if channel_energy < self._threshold + 0.1:
                    self._selectedchannels.append(channel)

    def calchannellabel(self, preci=1):
        i = 0
        for channel in self._selectedchannels:
            key1 = (self.voidid_label[channel.start].label,  self.voidid_label[channel.end].label,
                     round(channel.radii, preci))
            key2 = (self.voidid_label[channel.end].label,  self.voidid_label[channel.start].label,
                    round(channel.radii, preci))
            if key1 not in self.nonequalchannels.keys() and key2 not in self.nonequalchannels.keys():
                self.nonequalchannels[key1] = i
                i += 1
        for channel in self._selectedchannels:
            if (self.voidid_label[channel.start].label, self.voidid_label[channel.end].label,
                 round(channel.radii, preci)) in self.nonequalchannels.keys():
                key_nonequalchannel = (self.voidid_label[channel.start].label, self.voidid_label[channel.end].label,
                                        round(channel.radii, preci))
            elif (self.voidid_label[channel.end].label, self.voidid_label[channel.start].label,
                                        round(channel.radii, preci)) in self.nonequalchannels.keys():
                key_nonequalchannel = (self.voidid_label[channel.end].label, self.voidid_label[channel.start].label,
                                        round(channel.radii, preci))
            else:
                raise ValueError("非等价路径片段计算错误")
            channel.label = self.nonequalchannels[key_nonequalchannel]

    def to_net(self, filename_cif):
        with open(filename_cif.split(".")[0]+'_selectvoidbyBVSE_network.net', 'w') as f:
        # with open("./non_equivalent_paths_outputs/" + filename_cif.split("/")[-1].split(".")[0] + '_selectvoidbyBVSE_network.net', 'w') as f:
            f.write('Interstitial table:\n')
            for void in self._selectedvoids:
                f.write(str(void.id)+"\t"+str(void.label)+"\t"+str(void.coord[0])+" "
                        + str(void.coord[1]) + " "+str(void.coord[2])+"\t"+str(void.radii)+"\t"+str(void.energy)+"\n")
            f.write('Connection table:\n')
            for channel in self._selectedchannels:
                f.write(str(channel.start) + "\t " + str(channel.end) + "\t " + str(channel.phase[0]) + " "
                        + str(channel.phase[1]) + " " + str(channel.phase[2]) + "\t "
                        + str(channel.coord[0]) + " "
                        + str(channel.coord[1]) + " " + str(channel.coord[2]) + "\t "
                        + str(channel.radii) + "\t "
                        + str(channel.dist) + "\t "
                        + str(channel.label) + "\n")
