'''
__author__ = "Penghui Mi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Penghui Mi"
__email__ = "phmim@shu.edu.cn"
__date__ = "Mar 19, 2019"
'''
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
import networkx as nx
from monty.io import zopen
import ase.spacegroup as spg
from cavd.local_environment import CifParser_new
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure


class Void(object):
    """
    把间隙抽象为Void类，包含id、类别、坐标、半径、能量等属性
    """
    def __init__(self):
        self.id = None
        self.label = None
        self.coord = None
        self.radii = None
        self.energy = None


class Channel(object):
    """
    把通道抽象为Channel类，包含开始间隙点id、结束间隙点id、晶相、瓶颈分数坐标、瓶颈尺寸等属性
    """
    def __init__(self):
        self.start = None
        self.end = None
        self.phase = None
        self.coord = None
        self.radii = None
        self.dist = None
        self.label = None
        self.energy = None


def load_struc(filename_cif):
    """
    从cif文件中读取结构，调用的是pymatgen库中Structure类方法
    :param filename_cif: cif文件名
    :return: pymatgen定义的Structure对象
    """
    structure = Structure.from_file(filename_cif)
    return structure


def load_bvse(filename_bvse):
    """
    从numpy二进制文件中读取数据
    :param filename_bvse: 保存BVSE能量的文件
    :return: 三维数组
    """
    energy_landscape = np.load(filename_bvse)
    return energy_landscape


def load_voids_channels_from_file(filename_cavd):
    """
    从CAVD计算出的NET文件中读取间隙、导通信息
    :param filename_cavd: 要读取的文件名
    :return: 返回voids 和channels 两个字典
    """
    voids_dict = {}
    channels_dict = {}
    flag_p = 0
    flag_n = 0
    file = open(filename_cavd, 'r')
    for line in file.readlines():
        if 'Interstitial' in line:
            flag_p = 1
            flag_n = 0
            continue
        if 'Connection' in line:
            flag_p = 0
            flag_n = 1
            continue
        if flag_p == 1:
            line = line.split()
            if len(line) > 3:
                void = Void()
                void.id = int(line[0])
                void.label = int(line[1])
                void.coord = [np.float64(line[2]), np.float64(line[3]), np.float64(line[4])]
                void.radii = np.float64(line[5])
                voids_dict[void.id] = void
        if flag_n == 1:
            line = line.split()
            if len(line) > 4:
                channel = Channel()
                channel.start = int(line[0])
                channel.end = int(line[1])
                channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                channel.coord = [np.float64(line[5]), np.float64(line[6]), np.float64(line[7])]
                channel.radii = np.float64(line[8])
                channels_dict[(channel.start, channel.end)] = channel
    return voids_dict, channels_dict


class MergeCluster(object):
    def __init__(self, voids_dict, channels_dict, structure, bvse,filename_cif, clusterradii=0.5, radii_it=0.4, radii_bn=0.4):
        """
        根据间隙簇的距离来合并间隙簇
        :param clusterradii: 寻找间隙簇的距离判据
        :param radii_it:  用于筛选间隙的半径阈值
        :param radii_bn:  用于筛选通道的半径阈值
        """
        self._struc = structure
        self._energy = bvse
        self._voids = {}
        self._channels = {}
        self._mergedvoids = []
        self._mergedchannels = []
        self._clusterradii = clusterradii
        self._clusters = []
        self.init_voids_and_channels(voids_dict, channels_dict, radii_it, radii_bn)
        self.find_clusters()
        self.handle_voidsandchannels()
        self.cal_void_label(filename_cif)

    @property
    def mergedvoids(self):
        """
        合并后的间隙列表
        """
        return self._mergedvoids

    @property
    def mergedchannels(self):
        """
        合并后的通道片段列表
        """
        return self._mergedchannels

    def get_absolute_dis(self, p1, p2):
        """
        在不考虑周期性的情况下，计算两点之间的距离
        :param p1: 分数坐标，例如[0.5，0.5，0.5]
        :param p2: 分数坐标
        :return: 两点之间的欧式距离，类型为float
        """
        coord1 = np.array(self.fac2cart(p1))
        coord2 = np.array(self.fac2cart(p2))
        diff = np.subtract(coord1, coord2)
        return np.linalg.norm(diff)

    def get_period_dis(self, p1, p2):
        """
        在考虑周期性的情况下，计算两点之间的距离
        :param p1: 分数坐标，例如[0.5，0.5，0.5]
        :param p2: 分数坐标
        :return:  两点之间的距离，考虑周期性
        """
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def fac2cart(self, coord):
        """
        分数坐标转换成笛卡尔坐标
        """
        return np.dot(coord, self._struc.lattice.matrix)

    def cart2fac(self, coord):
        """
        笛卡尔坐标转换成分数坐标
        """
        return np.dot(coord, np.linalg.inv(self._struc.lattice.matrix))

    def cal_point_energy(self, point):
        # 计算能量场中一点的能量
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z), self._energy)
        for i in range(len(point)):
            if point[i] < 0.0:
                point[i] += 1
            if point[i] > 1.0:
                point[i] -= 1
        return interpolating_function(np.array(point))[0]

    def init_voids_and_channels(self, voids_dict, channel_dict, radii_it, radii_bn):
        for void_id, void in voids_dict.items():
            void.energy = self.cal_point_energy(void.coord)
        small_voids = [void_id for void_id, void in voids_dict.items() if void.radii < radii_it and void.energy < 5]
        self._voids = {void_id: void for void_id, void in voids_dict.items() if void.radii > radii_it}
        self._channels = {channel_id: channel for channel_id, channel in channel_dict.items()
                          if channel.start not in small_voids and channel.end not in small_voids}
        self._channels = {channel_id: channel for channel_id, channel in self._channels.items()
                          if channel.radii > radii_bn}

    def find_clusters(self):
        """
        在所有间隙中找到每一个簇
        """
        coords = []
        pair_voids = []
        voids_list = [void for void_id, void in self._voids.items()]
        for void in voids_list:
            coords.append(self.fac2cart(void.coord))
        coord_tree = cKDTree(coords)
        all_pair_voids = [i for i in coord_tree.query_pairs(r=self._clusterradii)]
        for i in all_pair_voids:
            pair_voids.append([voids_list[i[0]].id, voids_list[i[1]].id])
        if len(pair_voids) > 0:
            graph_clusters = nx.Graph()
            for e in pair_voids:
                graph_clusters.add_edge(e[0], e[1])
            queue_clusters = []
            for sc in nx.connected_component_subgraphs(graph_clusters):
                queue_clusters.append(list(sc.nodes))
            while queue_clusters:
                subv_in = queue_clusters.pop(0)
                subv_out = []
                temp_coord = [0.0, 0.0, 0.0]
                for subv in subv_in:
                    temp_coord[0] += self.fac2cart(self._voids[subv].coord)[0]
                    temp_coord[1] += self.fac2cart(self._voids[subv].coord)[1]
                    temp_coord[2] += self.fac2cart(self._voids[subv].coord)[2]
                centre_coord = self.cart2fac([temp_coord[0]/len(subv_in),
                                              temp_coord[1] / len(subv_in),
                                              temp_coord[2] / len(subv_in)])
                for i in range(len(subv_in)):
                    if self.get_period_dis(self._voids[subv_in[i]].coord, centre_coord) > self._clusterradii+0.2:
                        subv_out.append(subv_in[i])
                for subv in subv_out:
                    subv_in.remove(subv)
                self._clusters.append(subv_in)
                if len(subv_out) > 1:
                    pair_subvout = []
                    for i in range(len(subv_out)):
                        for j in range(i+1, len(subv_out)):
                            if self.get_period_dis(self._voids[subv_out[i]].coord,
                                                   self._voids[subv_out[j]].coord) < self._clusterradii:
                                pair_subvout.append([subv_out[i], subv_out[j]])
                    if len(pair_subvout) > 0:
                        graph_subvout = nx.Graph()
                        for e in pair_subvout:
                            graph_subvout.add_edge(e[0], e[1])
                        for sc in nx.connected_component_subgraphs(graph_subvout):
                            queue_clusters.append(list(sc.nodes))
        #print("间隙簇的个数为", len(self._clusters))
        #print(self._clusters)

    def handle_voidsandchannels(self):
        mignet = nx.Graph()
        for void_id, void in self._voids.items():
            mignet.add_node(void.id, label=void.label, coord=void.coord, radii=void.radii)
        for e_id, e in self._channels.items():  # 添加边
            if e.start < e.end:
                phase1 = e.phase
                phase2 = [0, 0, 0]
                for i in range(3):
                    if phase1[i] != 0:
                        phase2[i] = -1 * phase1[i]
                # phase1表示从id小的到大的phase,phase2表示从id大的到小的phase
                mignet.add_edge(e.start, e.end, phase1=phase1, phase2=phase2, coord1=e.coord, radii1=e.radii,
                                dist1=e.dist)
        if len(self._clusters) > 0:
            for i in range(len(self._clusters)):
                minenergy = 20000
                centervoid_id = self._clusters[i][0]
                for void_id in self._clusters[i]:
                    if minenergy > self._voids[void_id].energy:
                        centervoid_id = void_id
                        minenergy = self._voids[void_id].energy
                center_void = Void()
                center_void.id = centervoid_id
                center_void.label = self._voids[centervoid_id].label
                center_void.coord = self._voids[centervoid_id].coord
                center_void.radii = self._voids[centervoid_id].radii
                center_void.energy = self._voids[centervoid_id].energy
                tempedges = []
                nearvoids = []
                for nearvoid in list(mignet.adj[centervoid_id].keys()):
                    if nearvoid not in self._clusters[i]:
                        nearvoids.append(nearvoid)
                        if centervoid_id < nearvoid:
                            start = centervoid_id
                            end = nearvoid
                        else:
                            end = centervoid_id
                            start = nearvoid
                        tempedges.append({"from": start, "to": end,
                                          "phase1": mignet[centervoid_id][nearvoid]['phase1'],
                                          "phase2": mignet[centervoid_id][nearvoid]['phase2'],
                                          "coord1": mignet[centervoid_id][nearvoid]['coord1'],
                                          "radii1": mignet[centervoid_id][nearvoid]['radii1'],
                                          "dist1": mignet[centervoid_id][nearvoid]['dist1']})
                for id in self._clusters[i]:
                    if id != center_void.id:
                        for nearvoid in list(mignet.adj[id].keys()):
                            if nearvoid not in self._clusters[i] and nearvoid not in nearvoids:
                                if centervoid_id < nearvoid:
                                    start = centervoid_id
                                    end = nearvoid
                                else:
                                    end = centervoid_id
                                    start = nearvoid
                                if id < nearvoid:
                                    ph_cen_nearvoid = mignet[id][nearvoid]['phase1']
                                    ph_nearvoid_cen = mignet[id][nearvoid]['phase2']
                                else:
                                    ph_cen_nearvoid = mignet[id][nearvoid]['phase2']
                                    ph_nearvoid_cen = mignet[id][nearvoid]['phase1']
                                if centervoid_id<nearvoid:
                                    ph1 = ph_cen_nearvoid
                                    ph2 = ph_nearvoid_cen
                                else:
                                    ph2 = ph_cen_nearvoid
                                    ph1 = ph_nearvoid_cen
                                tempedges.append({"from": start, "to": end,
                                                  "phase1": ph1,
                                                  "phase2": ph2,
                                                  "coord1": mignet[id][nearvoid]['coord1'],
                                                  "radii1": mignet[id][nearvoid]['radii1'],
                                                  "dist1": mignet[id][nearvoid]['dist1']})
                for void in self._clusters[i]:
                    mignet.remove_node(void)
                mignet.add_node(center_void.id, label=center_void.label, coord=center_void.coord,
                                radii=center_void.radii)
                for e in tempedges:
                    mignet.add_edge(e['from'], e['to'],  phase1=e['phase1'], phase2=e['phase2'],
                                    coord1=e['coord1'], radii1=e['radii1'], dist1=e['dist1'])
        for nd in mignet.nodes():
            tempvoid = Void()
            tempvoid.id = nd
            tempvoid.label = mignet.nodes[nd]['label']
            tempvoid.coord = mignet.nodes[nd]['coord']
            tempvoid.radii = mignet.nodes[nd]['radii']
            self._mergedvoids.append(tempvoid)
        for edge in mignet.edges():
            tempchannel1 = Channel()
            tempchannel2 = Channel()
            tempchannel1.start = edge[0]
            tempchannel1.end = edge[1]
            tempchannel2.end = edge[0]
            tempchannel2.start = edge[1]
            if edge[0] < edge[1]:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase1"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase2"]
            else:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase2"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase1"]
            tempchannel1.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel2.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel1.radii = mignet[edge[0]][edge[1]]["radii1"]
            tempchannel2.radii = mignet[edge[0]][edge[1]]["radii1"]
            tempchannel1.dist = mignet[edge[0]][edge[1]]["dist1"]
            tempchannel2.dist = mignet[edge[0]][edge[1]]["dist1"]
            self._mergedchannels.append(tempchannel1)
            self._mergedchannels.append(tempchannel2)

    @staticmethod
    def tag_sites(sitesym, scaled_positions, symprec=1e-5):
        scaled = np.around(np.array(scaled_positions, ndmin=2), 8)
        scaled %= 1.0
        scaled %= 1.0
        np.set_printoptions(suppress=True)
        tags = -np.ones((len(scaled),), dtype=int)
        tagdis = 100 * np.ones((len(scaled),), dtype=float)
        rot, trans = spg.spacegroup.parse_sitesym(sitesym)
        siteskdTree = cKDTree(scaled)
        for i in range(len(scaled)):
            if tags[i] == -1:
                curpos = scaled[i]
                sympos = np.dot(rot, curpos) + trans
                sympos %= 1.0
                sympos %= 1.0
                sympos = np.unique(np.around(sympos, 8), axis=0)
                min_dis, min_ids = siteskdTree.query(sympos, k=1)
                select = min_dis < symprec
                select_ids = min_ids[select]
                tags[select_ids] = i
                tagdis[select_ids] = min_dis[select]
        return tags, tagdis

    def cal_void_label(self, filename_cif, symprec=0.1):
        with zopen(filename_cif, "rt") as f:
            input_string = f.read()
        parser = CifParser_new.from_string(input_string)
        sitesym = parser.get_sym_opt()

        voids_positions = []
        for void in self.mergedvoids:
            voids_positions.append(void.coord)
        tags, tagdis = self.tag_sites(sitesym, voids_positions, symprec)
        for i in range(len(tags)):
            self.mergedvoids[i].label = tags[i]

    def to_net(self, filename):
        """
        将合并之后的间隙点和通道片段保存到NET文件中
        :param filename: output will be written to a file
        """
        with open(filename.split(".")[0]+'_mergecluster.net', 'w') as f:
        # with open("./non_equivalent_paths_outputs/" + filename.split("/")[-1].split(".")[0] + '_mergecluster.net', 'w') as f:
            f.write('Interstitial table:\n')
            for void in self.mergedvoids:
                f.write(str(void.id)+"\t"+str(void.label)+"\t "
                        + str(void.coord[0]) + " " + str(void.coord[1]) + " "+str(void.coord[2]) + "\t "
                        + str(void.radii)
                        + "\n")
            f.write('Connection table:\n')
            for channel in self.mergedchannels:
                f.write(str(channel.start) + "\t " + str(channel.end) + "\t " + str(int(channel.phase[0])) + " "
                        + str(int(channel.phase[1])) + " " + str(int(channel.phase[2])) + "\t "
                        + str(channel.coord[0]) + " "
                        + str(channel.coord[1]) + " " + str(channel.coord[2]) + "\t "
                        + str(channel.radii) + "\t " + str(channel.dist) + "\n")

    def get_clusters(self):
        """
        Find all interstice clusters in a periodic unit cell
        :return: list
        """
        clusters = []
        for cluster in self._clusters:
            cluster_temp = []
            for void_id in cluster:
                cluster_temp.append({"void_id": void_id, "void_coord": self._voids[void_id].coord})
            clusters.append(cluster_temp)
        return clusters