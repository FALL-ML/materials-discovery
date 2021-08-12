'''
__author__ = "Penghui Mi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Penghui Mi"
__email__ = "phmim@shu.edu.cn"
__date__ = "Mar 19, 2019"
'''
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.sites import PeriodicSite
from cavd_bvse.mergecluster import Void, Channel
import numpy as np
import networkx as nx
from scipy.interpolate import RegularGridInterpolator
import numpy.linalg as la
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import operator
import copy


class MigrationPaths(object):
    def __init__(self, struc, energy, moveion="Li"):
        self.struc = struc
        self.energy = energy
        self.__voids = None
        self.__channels = None
        self.moveion = moveion
        self.Mignet = None
        self.nodles = []
        self.nodles_pair = None
        self.labels = []
        self.initpoints_path = None
        self.__paths = None
        self.__pathsenergy = None

    def setvoids(self, voids_list):
        self.__voids = voids_list

    def setchannels(self, channels_list):
        self.__channels = channels_list

    @property
    def paths(self):
        return self.__paths

    @property
    def pathsenergy(self):
        return self.__pathsenergy

    def get_dis(self, p1, p2):
        temp_site1 = PeriodicSite('Ar', p1, self.struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self.struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def bulid_migrationnet(self):
        self.Mignet = nx.Graph()
        for n in self.__voids:  # 添加点
            self.Mignet.add_node(n.id, label=n.label, coord=n.coord)

        for e in self.__channels:  # 添加边
            if e.start < e.end:
                l1 = e.phase
                l2 = [-1*i for i in e.phase]
            else:
                l1 = [-1*i for i in e.phase]
                l2 = e.phase
            self.Mignet.add_edge(e.start, e.end, label=e.label, phase1=l1, phase2=l2)

    def read_cif_nodles(self):
        self.labels = []
        a = SpacegroupAnalyzer(self.struc, symprec=0.1)
        symm_structure = a.get_symmetrized_structure()
        # nodles_pair[1,1]={'ids':[3,5],'dis':3.89604}
        for i, sites in enumerate(symm_structure.equivalent_sites):
            for site in sites:
                if site.specie.symbol == self.moveion:
                    nodle = {}
                    nodle['fracoord'] = site.frac_coords
                    nodle['gridcoord'] = np.array([site.frac_coords[0] * (self.energy.shape[0]-1),
                                                   site.frac_coords[1] * (self.energy.shape[1]-1),
                                                  site.frac_coords[2] * (self.energy.shape[2]-1)])
                    mini_dis = 100
                    mini_id = -1
                    for void in self.__voids:
                        dis = self.get_dis(nodle['fracoord'], void.coord)
                        if dis < mini_dis:
                            mini_id = void.id
                            mini_dis = dis
                            nodle['label'] = void.label
                            if mini_dis < 0.2:
                                break
                    nodle['id_CAVD'] = mini_id
                    if nodle['label'] not in self.labels:
                        self.labels.append(nodle['label'])
                    self.nodles.append(nodle)

        self.nodles_pair = {}
        min_dis = {}
        for i in range(len(self.nodles) - 1):
            for j in range(i + 1, len(self.nodles)):
                dis = self.get_dis(self.nodles[i]['fracoord'], self.nodles[j]['fracoord'])
                if (self.nodles[i]['label'], self.nodles[j]['label']) not in self.nodles_pair:
                    self.nodles_pair[self.nodles[i]['label'], self.nodles[j]['label']] = [
                        {'fracoords': [self.nodles[i]['fracoord'], self.nodles[j]['fracoord']],
                         'ids_CAVD': [self.nodles[i]['id_CAVD'], self.nodles[j]['id_CAVD']], 'dis': dis}]
                    min_dis[self.nodles[i]['label'], self.nodles[j]['label']] = dis
                else:
                    if dis < min_dis[self.nodles[i]['label'], self.nodles[j]['label']]:
                        min_dis[self.nodles[i]['label'], self.nodles[j]['label']] = dis
                    self.nodles_pair[self.nodles[i]['label'], self.nodles[j]['label']].append(
                        {'fracoords': [self.nodles[i]['fracoord'], self.nodles[j]['fracoord']],
                         'ids_CAVD': [self.nodles[i]['id_CAVD'], self.nodles[j]['id_CAVD']], 'dis': dis})
        #print(self.nodles_pair)

    def find_endpoints(self):
        self.initpoints_path = {}
        for key, value in self.nodles_pair.items():
            paths_dict = {}
            for v in value:
                paths_single = list(nx.all_simple_paths(self.Mignet,
                                               source=v['ids_CAVD'][0], target=v['ids_CAVD'][1], cutoff=5))
                # paths_single = [nx.shortest_path(self.Mignet,
                                              # source=v['ids_CAVD'][0], target=v['ids_CAVD'][1])]
                if paths_single:
                    for path in paths_single:
                        if len(path) < 6:
                            temppath = {}
                            pathkey = []
                            path_phase = [0,0,0]
                            for i in range(len(path)-1):
                                # 计算路径片段的label
                                pathkey.append(self.Mignet[path[i]][path[i+1]]['label'])
                                if path[i] < path[i+1]:
                                    path_phase[0] += self.Mignet[path[i]][path[i+1]]['phase1'][0]
                                    path_phase[1] += self.Mignet[path[i]][path[i + 1]]['phase1'][1]
                                    path_phase[2] += self.Mignet[path[i]][path[i + 1]]['phase1'][2]
                                else:
                                    path_phase[0] += self.Mignet[path[i]][path[i + 1]]['phase2'][0]
                                    path_phase[1] += self.Mignet[path[i]][path[i + 1]]['phase2'][1]
                                    path_phase[2] += self.Mignet[path[i]][path[i + 1]]['phase2'][2]

                            a = True  # a用来判断是否含义晶格位点
                            path_voidlabel = [self.Mignet.node[path[0]]["label"]]
                            for i in range(1, len(path)-1):
                                lab = self.Mignet.node[path[i]]["label"]
                                if lab not in self.labels:
                                    path_voidlabel.append(lab)
                                else:
                                    a = False
                                    break
                            path_voidlabel.append(self.Mignet.node[path[-1]]["label"])
                            temppath['path_voidlabel'] = path_voidlabel
                            temppath['path_voidid'] = path
                            temppath['path_phase'] = path_phase

                            if a and tuple(pathkey) not in paths_dict.keys() and\
                                    tuple(pathkey[::-1]) not in paths_dict.keys():
                                paths_dict[tuple(pathkey)] = temppath

                            if a and tuple(pathkey) in paths_dict.keys():
                                if operator.eq(temppath['path_phase'], [0, 0, 0]) and \
                                        not (operator.eq(paths_dict[tuple(pathkey)]['path_phase'], [0, 0, 0])):
                                    paths_dict[tuple(pathkey)] = temppath

                            if a and tuple(pathkey[::-1]) in paths_dict.keys():
                                if operator.eq(temppath['path_phase'], [0, 0, 0]) and \
                                        not (operator.eq(paths_dict[tuple(pathkey[::-1])]['path_phase'], [0, 0, 0])):
                                    paths_dict[tuple(pathkey[::-1])] = temppath

            for key, value in paths_dict.items():
                if value['path_voidid'][0] < value['path_voidid'][-1]:
                    self.initpoints_path[(value['path_voidid'][0], value['path_voidid'][-1])] ={
                        "path_voidid": value['path_voidid'],"phase" :value["path_phase"] }
                else:
                    self.initpoints_path[(value['path_voidid'][-1], value['path_voidid'][0])] = {
                        "path_voidid": value['path_voidid'][::-1], "phase" :[-1*i for i in value["path_phase"]] }


    def cal_equal_path(self):
        self.__paths = []
        for initpoints in self.initpoints_path.values():
            if operator.eq(initpoints["phase"], [0, 0, 0]):
                p = self.calincellpath([initpoints['path_voidid'][0], initpoints['path_voidid'][-1]], count_images=21)

            else:
                temp_inipoints = [initpoints['path_voidid'][0], initpoints['path_voidid'][-1]] + initpoints["phase"]
                p = self.caloutcellpath(temp_inipoints, count_images=21)
            if p is not None:
                self.__paths.append(p)
            #print(p)

    def calpointenergy(self, point):
        # 计算能量场中一点的能量
        x = np.linspace(0, self.energy.shape[0]-1, self.energy.shape[0])
        y = np.linspace(0, self.energy.shape[1]-1, self.energy.shape[1])
        z = np.linspace(0, self.energy.shape[2]-1, self.energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z), self.energy)
        return interpolating_function(np.array(point))

    def pathenergy(self, p):
        self.energy = self.energy - np.amin(self.energy)
        energy_path = []
        for point_of_p in p:
            point_temp = [0.0, 0.0, 0.0]
            for i in range(len(point_temp)):
                if point_of_p[i] < 0:
                    point_temp[i] = (point_of_p[i]+1)*(self.energy.shape[i]-1)
                elif point_of_p[i] < 1:
                    point_temp[i] = point_of_p[i] * (self.energy.shape[i]-1)
                else:
                    point_temp[i] = (point_of_p[i]-1) * (self.energy.shape[i] - 1)
            energy_point_temp = self.calpointenergy(point_temp)
            energy_path.append(energy_point_temp)
        dis_path = [0.0]
        tol = 0
        for i in range(len(p)-1):
            dist = self.get_dis(p[i], p[i+1])
            tol += dist
            dis_path.append(tol)
        migs = np.zeros(shape=(len(dis_path), 2))
        for i in range(len(energy_path)):
            migs[i][1] = energy_path[i]
            migs[i][0] = dis_path[i]
        return migs

    def calpathsenergy(self):
        # 输出分别为横坐标 ，纵坐标
        vaipaths = self.__paths
        self.__pathsenergy = []
        for edge in vaipaths:
            eneg = self.pathenergy(edge)
            self.__pathsenergy.append(eneg)

    def showenergy(self, filename_cif):
        for i in range(len(self.__pathsenergy)):
            xcoords = []
            ycoords = []
            for j in range(len(self.__pathsenergy[i])):
                xcoords.append(self.__pathsenergy[i][j][0])
                ycoords.append(self.__pathsenergy[i][j][1])
            poly = np.polyfit(xcoords, ycoords, deg=12)  # 最小二乘法多项式拟合
            z = np.polyval(poly, xcoords)
            plt.figure(figsize=(6, 4))  # 创建绘图对象
            plt.plot(xcoords, z, linewidth=3)
            plt.scatter(xcoords, ycoords, color='k', marker='o')
            plt.xlabel("Reaction Coordinate ")  # X轴标签
            plt.ylabel("Energy")  # Y轴标签
            path = filename_cif.split(".")[0] + str(i) + '.png'
            plt.savefig(path)
            plt.show()  # 显示图

    @staticmethod
    def pathfind(start, end, V, n_images=21, dr=None, h=0.1, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6):
        # Set parameters
        if not dr:
            dr = np.array([1.0 / (V.shape[0]-1), 1.0 / (V.shape[1]-1), 1.0 / (V.shape[2]-1)])
        else:
            dr = np.array(dr, dtype=float)
        keff = k * dr * n_images
        h0 = h
        # 初使化 string
        g1 = np.linspace(0, 1, n_images)  # 创建等差数列
        s0 = start  # 初使结构的迁移离子的坐标
        s1 = end  # 最终结构的迁移离子的坐标
        s = np.array([g * (s1 - s0) for g in g1]) + s0  # s是一个3*n的矩阵
        ds = s - np.roll(s, 1, axis=0)  # np.roll函数的意思是将s，沿着axis的方向，滚动1个长度，相当于把最后一个元素放到第一个
        ds[0] = (ds[0] - ds[0])  # 把第一个元素置为（0，0，0）
        ls = np.cumsum(la.norm(ds, axis=1))  # norm求ds的范数  cumsum计算轴向元素累加和，返回由中间结果组成的数组
        ls = ls / ls[-1]  # 归一化
        fi = interp1d(ls, s, axis=0)  # 插值
        s = fi(g1)
        #print(s)
        # Evaluate initial distances (for elastic equilibrium)
        ds0_plus = s - np.roll(s, 1, axis=0)  # 正向
        ds0_minus = s - np.roll(s, -1, axis=0)  # 负向
        ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
        ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])
        dV = np.gradient(V)  # 计算梯度

        # Evolve string
        for step in range(0, max_iter):
            #if step > min_iter:
                #h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)  # 逐步衰减步长以减少震荡
            #else:
            h = h0
            # Calculate forces acting on string
            d = V.shape
            s0 = s
            edV = np.array([[dV[0][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[0],

                             dV[1][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[1],

                             dV[2][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[2]] for pt in s])

            # Update according to force due to potential and string elasticity
            ds_plus = s - np.roll(s, 1, axis=0)
            ds_minus = s - np.roll(s, -1, axis=0)
            ds_plus[0] = (ds_plus[0] - ds_plus[0])
            ds_minus[-1] = (ds_minus[-1] - ds_minus[-1])

            Fpot = edV
            Fel = keff * (la.norm(ds_plus) - la.norm(ds0_plus)) * (ds_plus / la.norm(ds_plus))
            Fel += keff * (la.norm(ds_minus) - la.norm(ds0_minus)) * (ds_minus / la.norm(ds_minus))
            s = s - h * (Fpot + Fel)
            # 每次迭代保持两端值保持不变
            s[0] = s0[0]
            s[-1] = s0[-1]
            # Reparametrize string            #更新参数
            ds = s - np.roll(s, 1, axis=0)
            ds[0] = (ds[0] - ds[0])
            ls = np.cumsum(la.norm(ds, axis=1))
            ls = ls / ls[-1]
            fi = interp1d(ls, s, axis=0)
            s = fi(g1)
            tol = la.norm((s - s0) * dr) / n_images / h
            if (tol > 1e9):
                s=[[0,0,0]]
                break
            if (step > min_iter and tol < max_tol):
                # print("Converged at step {}".format(step))
                break
            # if (step % 100 == 0):
            # print ("Step {} - ds = {}".format(step, tol))
        return s

    def calincellpath(self, initpoints, count_images=9):
        # initpoints是一个列表，initial[0]表示起点坐标，initial[1]表示终点坐标，initial[2，3，4]表示终点坐标的晶相
        start = np.array([0.0, 0.0, 0.0])
        end = np.array([0.0, 0.0, 0.0])
        for void in self.__voids:
            if initpoints[0] == void.id:
                start[0] = void.coord[0] * (self.energy.shape[0] - 1)
                start[1] = void.coord[1] * (self.energy.shape[1] - 1)
                start[2] = void.coord[2] * (self.energy.shape[2] - 1)
            if initpoints[1] == void.id:
                end[0] = void.coord[0] * (self.energy.shape[0] - 1)
                end[1] = void.coord[1] * (self.energy.shape[1] - 1)
                end[2] = void.coord[2] * (self.energy.shape[2] - 1)

        p = self.pathfind(start, end, self.energy, n_images=count_images,
                              dr=[self.struc.lattice.a / (self.energy.shape[0]-1),
                                  self.struc.lattice.b / (self.energy.shape[1]-1),
                                  self.struc.lattice.c / (self.energy.shape[2]-1)], h=0.02, k=0.17, min_iter=100,
                              max_iter=10000, max_tol=5e-6)

        if len(p) > 1:
                istrue = 0
                for p1 in p:
                    p1[0] = round(p1[0] / (self.energy.shape[0]-1), 4)
                    p1[1] = round(p1[1] / (self.energy.shape[1]-1), 4)
                    p1[2] = round(p1[2] / (self.energy.shape[2]-1), 4)
                    if p1[0] < -1.1 or p1[0] > 2 or p1[1] < -1.1 or p1[1] > 2 or p1[2] < -1.1 or p1[2] > 2:
                        istrue = 1
                        break
                if istrue == 0:
                    return p

    def caloutcellpath(self, initpoints, count_images=11):
        a2 = np.concatenate((self.energy, self.energy, self.energy), axis=0)
        a3 = np.concatenate((a2, a2, a2), axis=1)
        expan_arr = np.concatenate((a3, a3, a3), axis=2)
        start = np.array([0.0, 0.0, 0.0])
        end = np.array([0.0, 0.0, 0.0])
        for void in self.__voids:
            if initpoints[0] == void.id:
                start[0] = void.coord[0] * (self.energy.shape[0] - 1) + self.energy.shape[0]
                start[1] = void.coord[1] * (self.energy.shape[1] - 1) + self.energy.shape[1]
                start[2] = void.coord[2] * (self.energy.shape[2] - 1) + self.energy.shape[2]
            if initpoints[1] == void.id:
                end[0] = void.coord[0] * (self.energy.shape[0] - 1) +(1+initpoints[2])*(self.energy.shape[0])
                end[1] = void.coord[1] * (self.energy.shape[1] - 1) +(1+initpoints[3])*(self.energy.shape[1])
                end[2] = void.coord[2] * (self.energy.shape[2] - 1) +(1+initpoints[4])*(self.energy.shape[2])
        p = self.pathfind(start, end, expan_arr, n_images=count_images,
                          dr=[self.struc.lattice.a / (self.energy.shape[0]-1),
                              self.struc.lattice.b / (self.energy.shape[1]-1),
                              self.struc.lattice.c / (self.energy.shape[2]-1)], h=0.01, k=0.17, min_iter=100,
                          max_iter=10000, max_tol=5e-6)
        if len(p) > 1:
            istrue = 0
            for p1 in p:
                p1[0] = round((p1[0]-self.energy.shape[0]+1) / self.energy.shape[0], 4)
                p1[1] = round((p1[1]-self.energy.shape[1]+1) / self.energy.shape[1], 4)
                p1[2] = round((p1[2]-self.energy.shape[2]+1) / self.energy.shape[2], 4)
                if p1[0] < -1.1 or p1[0] > 2 or p1[1] < -1.1 or p1[1] > 2 or p1[2] < -1.1 or p1[2] > 2:
                    istrue = 1
                    break
            if istrue == 0:
                return p

    def savedata(self, filename_cif):
        with open(filename_cif.split(".")[0] + '_nonequalpaths.txt', 'w') as f:
        # with open("./non_equivalent_paths_outputs/" + filename_cif.split("/")[-1].split(".")[0] + '_nonequalpaths.txt', 'w') as f:
            for i in range(len(self.__paths)):
                f.write('nonequalpath id:' + str(i) + '\n')
                for j in range(len(self.__paths[i])):
                    f.write(str(j)+"\t"+str(self.__paths[i][j][0])+'  ' + str(self.__paths[i][j][1])+'  '
                            + str(self.__paths[i][j][2]) + '\t' + str(self.__pathsenergy[i][j][0])+'  '
                            + str(self.__pathsenergy[i][j][1])+'\n ')

