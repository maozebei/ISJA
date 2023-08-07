# -*- coding: utf-8 -*-
"""
Created on 2022/7/13 10:11

@author: MZB
"""

import numpy as np
import os
import time
import xml.etree.ElementTree as ET
import xml.dom.minidom as MD
import matplotlib.pyplot as plt
import matplotlib.tri as tri
# 工作目录
os.chdir(r'./temp')


## -------------------------------------------------- ##
input_name      = 'Job-demo'                # abaqus完整的inp文件名
study_name      = 'sdy_buckle'              # 当前任务名称
initial_name    = 'study'                   # moldflow空任务文件名(序列改为填充+保压)
stressfile      = 'stress_tensor.txt'       # 应力张量

# 注塑口控制
inj_type        = 3                         # 1:节点号  2:绝对坐标  3:归一化坐标
inj_method      = 0                         # 0:手动定义  1:随机生成
num_group       = [0]

# 生成注塑口位置参数
size_x          = 40                        # x轴分块
size_y          = 40                        # y轴分块
num_inj         = 1                         # 注塑口数量

# 材料及线程
mat_ID          = 13970                     # 35%CF+PPS
thread_MF       = 1                         # MoldFlow并行线程数
thread_ABQ      = 1                         # Abaqus并行线程数

if_C_W          = 0                         # 是否画图  0:不画图  1:无应力加权  2:应力加权

# /parameter
## ================================================== ##



## -------------------------------------------------- ##
# main

def main():

    # 前处理
    part, Others, mesh_type, NoE = Pre_processing(input_name, initial_name, study_name)

    # 注塑口参数化
    inj_group = Injection_Parametric(size_x, size_y, num_inj, num_group)
    (inj_location_group, inj_vector_group) = inj_group

    # 主循环
    for group in range(len(inj_location_group)):

        # 模流分析
        T1 = time.time()
        inj_loc = inj_location_group[group]
        inj_vec = inj_vector_group[group]
        Moldflow_Analysis(group, study_name, initial_name, mesh_type, 
                          inj_type, inj_loc, inj_vec, mat_ID, thread_MF)
        T2 = time.time()
        print('moldflow time: %f s' % (T2-T1))
        # 结构分析
        Abaqus_Analysis(group, study_name, part, NoE, Others, thread_ABQ)
    
    T3 = time.time()
    print('abaqus time: %f s' % (T3-T2))
    # 后处理
    Post_processing(num_group, NoE, stressfile, if_C_W)
    T4 = time.time()
    print('post time: %f s' % (T4-T3))

# /main
## ================================================== ##



## -------------------------------------------------- ##
# 分析模块

# 前处理
def Pre_processing(input_name, initial_name, study_name):

    # 读取inp文件信息
    Parts, Others = ReadInputFile(input_name)

    # 检查part是否只包含一个part, 输出单元坐标文件
    part, mesh_type = CheckInput(Parts)

    # 单元数量
    NoE = len(part.Element)

    # 转为bdf格式
    bdf_name = initial_name + '.bdf'
    WriteBDFfile(part, bdf_name)

    # 若存在屈曲特征值储存文件，删除
    try: os.remove('buckle_eigenvalue.txt')
    except: pass

    # 删除前一步计算留下的结果文件
    name = os.listdir(os.getcwd())
    for nm in name:
        if study_name in nm:
            os.remove(nm)

    return part, Others, mesh_type, NoE


# 注塑口位置随机参数化模块
def Injection_Parametric(size_x, size_y, num_inj, num_group):

    # 将设计域的点划分为 size_x+1行 size_y+1列
    P_num = (size_x + 1) * (size_y + 1)

    # 随机挑选注塑点
    P_ID = list()
    for inj in num_group:
        P_ID.append([inj])

    # 将P_ID转为归一化坐标
    inj_location = list()
    inj_vector = list()
    for i in range(len(P_ID)):
        inj_location_group = list()
        inj_vector_group = list()
        for j in range(len(P_ID[i])):
            x = P_ID[i][j] % (size_x + 1) * (1.0 / size_x)
            y = P_ID[i][j] // (size_x + 1) * (1.0 / size_y)
            z = 0.0
            inj_location_group.append([x, y, z])
            inj_vector_group.append([0.0, 0.0, 1.0])
        inj_location.append(inj_location_group)
        inj_vector.append(inj_vector_group)
            
    # 以[0 1]形式输出注塑口向量
    inj_total = list()
    for i in range(len(P_ID)):
        P_v = np.zeros(P_num, dtype=np.int32)
        for j in range(len(P_ID[i])):
            P_v[P_ID[i][j]] = 1
        inj_total.append(P_v)
    
    # 以[0 1]形式将注塑口矩阵存入txt
    np.savetxt('injection_location.txt', np.array(inj_total), delimiter=' ', fmt='%d')

    return (inj_location, inj_vector)


# 模流分析
def Moldflow_Analysis(group, study_name, initial_name, mesh_type, 
                      inj_type, inj_location, inj_vector, mat_ID, thread_MF):

    # 创建xml文件
    bdf_name = initial_name + '.bdf'
    CreatXML(mesh_type, bdf_name, inj_type, inj_location, inj_vector, mat_ID, thread_MF)

    # 通过studymod程序更改sdy
    os.system('studymod %s.sdy %s %s.xml' % (initial_name, study_name, initial_name))

    # 注塑模拟计算
    sdy = study_name
    if mesh_type == 1: os.system('flow -fill -pack -fiber -stress -runb -output %s~3 %s.sdy' % (sdy, sdy))
    if mesh_type == 2: os.system('mhb3d -fill -pack -output %s~3 %s.sdy' % (sdy, sdy))

    # 提取纤维取向张量
    if mesh_type == 1:
        os.system('studyrlt %s.sdy -xml 4000' % sdy)
        RenameFile('%s.xml' % sdy, 'group%d-fiber.xml' % group)
        os.system('studyrlt %s.sdy -xml 4020' % sdy)
        RenameFile('%s.xml' % sdy, 'group%d-E1.xml' % group)
        os.system('studyrlt %s.sdy -xml 4030' % sdy)
        RenameFile('%s.xml' % sdy, 'group%d-E2.xml' % group)
        os.system('studyrlt %s.sdy -xml 4040' % sdy)
        RenameFile('%s.xml' % sdy, 'group%d-G12.xml' % group)
        os.system('studyrlt %s.sdy -xml 4050' % sdy)
        RenameFile('%s.xml' % sdy, 'group%d-Nu.xml' % group)
    # if mesh_type == 2: os.system('studyrlt %s -xml 4009' % study_name)
    
    return


# 结构分析
def Abaqus_Analysis(group, study_name, part, NoE, Others, thread_ABQ):

    # 读取纤维取向结果和力学常数，输出inp附属文件 和 计算用的inp文件
    newinp_name = study_name + '.inp'
    InputIncludeFile(group, part, NoE)
    NewInputFile(part, Others, newinp_name)

    # abaqus计算
    os.system('abaqus job=%s cpus=%d int' % (newinp_name, thread_ABQ))

    # abaqus提取一阶屈曲特征值
    os.system('abaqus cae nogui=abaquspost_buckle.py')
    return


# 后处理
def Post_processing(num_group, NoE, stressfile, if_C_W):

    f = open('similar.txt', 'w')
    f2 = open('fiberorien_tensor.txt', 'w')

    for group in num_group:

        fiberfile = 'group%d-fiber.xml' % group
        coordfile = 'element_coordinate.txt'
        figname = 'similar.tif'

        # 计算单元相似度
        Ele_similar, Ele_similar_w, Fiber = SimilarCalculate(NoE, fiberfile, stressfile, 1)
        # 读取单元
        Ele_coord = np.loadtxt(coordfile)

        # 输出平均相似度和纤维取向张量
        f.write('%f %f\n' % (np.mean(Ele_similar), np.mean(Ele_similar_w)))
        for fb in Fiber:
            f2.write('%f %f %f %f %f %f ' % tuple(fb))
        f2.write('\n')

        # # 画相似度云图
        # if if_C_W == 1: plot_contour(Ele_coord, Ele_similar, figname)
        # elif if_C_W == 2: plot_contour(Ele_coord, Ele_similar_w, figname)
    
    f.close()
    f2.close()
    
    return

# /分析模块
## ================================================== ##




## -------------------------------------------------- ##
# Class 

# Part结构体
class iParts(object):
    def __init__(self, partname):
        self.Name = partname
        self.Node = list()
        self.Element = dict()
        self.ElementType = dict()
        self.NodeSet = dict()
        self.ElementSet = dict()
        self.Section = list()
        self.ifNodeX = list()
    
    def NodeRelabel(self):
        # 节点标号重新排列
        count = 0
        NodeID = dict()
        temp = list()
        # 更新节点，标号对应关系存在NodeID
        for i in range(len(self.Node)):
            if not self.ifNodeX[i]:
                count += 1
                NodeID[i+1] = count
                temp.append(self.Node[i])
            else:
                NodeID[i+1] = None
        self.Node = temp
        self.ifNodeX = [False] * len(temp)
        # 更新NodeSet中的节点号
        for s in self.NodeSet:
            temp = list()
            for n in self.NodeSet[s]:
                if NodeID[n]: temp.append(NodeID[n])
            self.NodeSet[s] = temp
        # 更新Element中的节点号
        for el in self.Element:
            temp = list()
            for n in self.Element[el]:
                if NodeID[n]: temp.append(NodeID[n])
            self.Element[el] = temp

    def ElementRelabel(self):
        # 单元号重新排列
        count = 0
        ElmtID = dict()
        temp1 = dict()
        temp2 = dict()
        # 更新单元号，标号对应关系存在ElmtID
        for el in self.Element:
            count += 1
            ElmtID[el] = count
            temp1[count] = self.Element[el]
            temp2[count] = self.ElementType[el]
        self.Element = temp1
        self.ElementType = temp2
        # 更新ElementSet中的单元号
        for s in self.ElementSet:
            temp = list()
            for el in self.ElementSet[s]:
                temp.append(ElmtID[el])
            self.ElementSet[s] = temp

# Material结构体
class iMaterials(object):
    def __init__(self, materialname):
        self.Name = materialname
        self.Elastic = np.zeros(0)
        self.Density = 0.
        self.Plastic = np.zeros(0)

# 分割单元结构体
class Element_Sp(object):
    def __init__(self, ElementType):
        self.ElementType = ElementType
        self.NodeID = list()
        self.NodeCoord = list()
        self.ifNodeX = list()
        self.Edge = list()      # [{n1_ID:Coord, n2_ID:Coord}, ...]
        self.NewNode = list()
        self.NoDeleteNode = 0   # 应删除的节点数
        self.NoNewNode = 0      # 应新增的节点数

# Finish Class
## ================================================== ##



## -------------------------------------------------- ##
# 读取inp文件信息

# 读取Input文件信息
def ReadInputFile(inp_name):
    # 从文件中读取信息
    filename = inp_name + '.inp'
    f = open(filename, 'r')
    InpLines = f.readlines()
    f.close()
    # 循环 读行
    NoL = 0
    Parts = dict()
    Materials = dict()
    Others = list()
    while NoL < len(InpLines):
        # 若字符串中有**则忽略(注释), 要求文件最后一行为空行
        lines = InpLines[NoL]
        line = lines[:lines.find('**')].strip().split(', ')
        # 根据抬头选择处理数据的函数
        if line[0] == '*Part':
            part, NoL = InputFilter_Part(InpLines, NoL)
            Parts[part.Name] = part
        elif line[0] == '*Material':
            material, NoL = InputFilter_Material(InpLines, NoL)
            Materials[material.Name] = material
        else:
            Others = InputFilter_Others(lines, Others)
            NoL += 1
    return Parts, Others

# 储存Part信息
def InputFilter_Part(InpLines, NoL):
    part = iParts(InpLines[NoL].strip().split(', ')[1][5:])
    while '*End Part' not in InpLines[NoL]:
        if NoL >= len(InpLines):
            break
        lines = InpLines[NoL]
        NoL += 1
        line = lines[:lines.find('**')].strip().split(', ')
        if line[0] == '':
            continue
        elif line[0] == '*Node':
            node, NoL = InputFilter_Part_Node(InpLines, NoL)
            part.Node = node
        elif line[0] == '*Element':
            element, NoL = InputFilter_Part_Element(InpLines, NoL)
            for el in element:
                part.ElementType[el[0]] = el[1]
                part.Element[el[0]] = el[2:]
        elif line[0] == '*Nset':
            nset, node, NoL = InputFilter_Part_Nset(InpLines, NoL)
            part.NodeSet[nset] = node
        elif line[0] == '*Elset':
            elset, element, NoL = InputFilter_Part_Elset(InpLines, NoL)
            part.ElementSet[elset] = element
        elif '*' in line[0] and 'Section' in line[0]:
            section, NoL = InputFilter_Part_Section(InpLines, NoL)
            part.Section.append(section)
    return part, NoL+1

# 读取part中的node
def InputFilter_Part_Node(InpLines, NoL):
    node = list()
    while '*' not in InpLines[NoL]:
        if NoL >= len(InpLines):
            break
        line = InpLines[NoL].strip().split(', ')
        temp = list()
        for l in line[1:]:
            temp.append(float(l))
        node.append(temp)
        NoL += 1
    return node, NoL

# 读取part中的element
def InputFilter_Part_Element(InpLines, NoL):
    element = list()
    type = InpLines[NoL-1].strip().split(', ')[1][5:]
    while '*' not in InpLines[NoL]:
        if NoL >= len(InpLines):
            break
        line = InpLines[NoL].strip().split(', ')
        temp = [int(line[0]), type]
        for l in line[1:]:
            temp.append(int(l))
        element.append(temp)
        NoL += 1
    return element, NoL

# 读取part中的节点set
def InputFilter_Part_Nset(InpLines, NoL):
    line0 = InpLines[NoL-1].strip().split(', ')
    nset = line0[1][5:]
    node = list()
    while '*' not in InpLines[NoL]:
        if NoL >= len(InpLines):
            break
        line = InpLines[NoL].strip().split(', ')
        if len(line) == 1 and line[0][-1] == ',':
            line[0] = line[0][:-1]
        if line0[-1] == 'generate':
            node += list(range(int(line[0]), int(line[1])+1, int(line[2])))
        else:
            for l in line:
                node.append(int(l))
        NoL += 1
    return nset, node, NoL

# 读取part中的单元set
def InputFilter_Part_Elset(InpLines, NoL):
    line0 = InpLines[NoL-1].strip().split(', ')
    elset = line0[1][6:]
    element = list()
    while '*' not in InpLines[NoL]:
        if NoL >= len(InpLines):
            break
        line = InpLines[NoL].strip().split(', ')
        if len(line) == 1 and line[0][-1] == ',':
            line[0] = line[0][:-1]
        if line0[-1] == 'generate':
            element += list(range(int(line[0]), int(line[1])+1, int(line[2])))
        else:
            for l in line:
                element.append(int(l))
        NoL += 1
    return elset, element, NoL

# 读取part中的Section
def InputFilter_Part_Section(InpLines, NoL):
    line0 = InpLines[NoL-1].strip().split(', ')
    stype = line0[0][1:-8]
    elset = line0[1][6:]
    material = line0[2][9:]
    section = [stype, elset, material]
    if stype == 'Shell':
        line = InpLines[NoL].strip().split(', ')
        section += [float(line[0]), int(line[1])]
    return section, NoL+1

# 储存Material信息
def InputFilter_Material(InpLines, NoL):
    material = iMaterials(InpLines[NoL].strip().split(', ')[1][5:])
    while NoL < len(InpLines)-1:
        NoL += 1
        lines = InpLines[NoL]
        line = lines[:lines.find('**')].strip().split(', ')
        # Elastic
        if line[0] == '*Elastic':
            NoL += 1
            lines = InpLines[NoL]
            line = lines[:lines.find('**')].strip().split(', ')
            temp = [float(line[0]), float(line[1])]
            material.Elastic = np.array(temp)
        # Density
        elif line[0] == '*Density':
            NoL += 1
            lines = InpLines[NoL]
            line = lines[:lines.find('**')].strip().split(',')
            material.Density = float(line[0])
        # Plastic
        elif line[0] == '*Plastic':
            NoL += 1
            temp = list()
            while '*' not in InpLines[NoL]:
                line = InpLines[NoL].strip().split(', ')
                temp.append([float(line[0]), float(line[1])])
                NoL += 1
                if NoL > len(InpLines)-1:
                    break
            material.Plastic = np.array(temp)
        else:
            break
    return material, NoL+1

# 储存其他信息
def InputFilter_Others(line, Others):
    if line[:2] == '**': return Others
    if '*Heading' in line: return Others
    if '*Preprint' in line: return Others
    Others.append(line)
    return Others

# 检查INP文件中的part数量, 输出单元位置信息
def CheckInput(Parts):
    if len(Parts) > 1:
        print('Error: input文件中只能包含1个part')
    elif len(Parts) < 1:
        print('Error: input文件未找到part')
    else:
        partname = list(Parts.keys())[0]
    if Parts[partname].Section[0][0] == 'Shell': mesh_type = 1
    elif Parts[partname].Section[0][0] == 'Solid': mesh_type = 2
    part = Parts[partname]
    Elcoord = np.zeros([len(part.Element), 3])
    for i in part.Element:
        nodes = part.Element[i]
        coord = np.zeros(3)
        count = 0
        for n in nodes:
            coord += np.array(part.Node[n-1])
            count += 1
        Elcoord[i-1] += coord/count
    np.savetxt('element_coordinate.txt', Elcoord, delimiter=' ')
    return part, mesh_type

# /读取inp文件信息
## ================================================== ##



## -------------------------------------------------- ##
# abaqus输出文件转化为moldflow计算文件

# 输出nastran bdf格式文件
def WriteBDFfile(Part, filename):
    
    f = open(filename,'w')

    for i in range(len(Part.Node)):
        f.write('GRID*    %15d                 ' % (i+1))
        f.write('%.9e %.9e\n' % (Part.Node[i][0],Part.Node[i][1]))
        f.write('*        %.9e\n' % Part.Node[i][2])

    NoS = 0
    for section in Part.Section:
        NoS += 1
        # 壳体
        if section[0] == 'Shell':
            f.write('PSHELL*  %15d %15d ' % (NoS, 1))
            f.write('%.9e %15d\n' % (section[3], 1))
            f.write('*                                      1\n')
            for el in Part.ElementSet[section[1]]:
                nodes = Part.Element[el]
                # 三角形单元
                if Part.ElementType[el] == 'S3':
                    f.write('CTRIA3*  %15d %15d %15d %15d\n' % (el, NoS, nodes[0], nodes[1]))
                    f.write('*        %15d\n' % nodes[2])
                # 四边形单元
                elif Part.ElementType[el] == 'S4R':
                    f.write('CQUAD4*  %15d %15d %15d %15d\n' % (el, NoS, nodes[0], nodes[1]))
                    f.write('*        %15d %15d\n' % (nodes[2], nodes[3]))
                else:
                    print('Error: 壳体单元类型应当为[一阶三角形单元(S3) or 一阶四边形单元(S4R)]')
        # 实体
        if section[0] == 'Solid':
            f.write('PSOLID*  %15d %15d\n' % (NoS, 0))
            f.write('*\n')
            for el in Part.ElementSet[section[1]]:
                nodes = Part.Element[el]
                # 四面体单元
                if Part.ElementType[el] == 'C3D4':
                    f.write('CTETRA*  %15d %15d %15d %15d\n' % (el, NoS, nodes[0], nodes[1]))
                    f.write('*        %15d %15d\n' % (nodes[2], nodes[3]))
                else:
                    print('Error: 实体单元类型应当为[一阶四面体单元(C3D4)]')

    f.close()
    return

# 创建xml文件
def CreatXML(meshtype, bdf_name, inj_type, inj_location, inj_vector, mat_ID, thread):
    # 创建dom文档
    IMPL = MD.getDOMImplementation()
    dom = IMPL.createDocument(None, 'StudyMod', None)

    # 创建根节点
    root = dom.documentElement
    root.setAttribute('title', 'Autodesk StudyMod')
    root.setAttribute('ver', '1.00')

    # 单位
    UnitSystem = dom.createElement('UnitSystem')
    text = dom.createTextNode('Metric')
    UnitSystem.appendChild(text)
    root.appendChild(UnitSystem)
    
    # 网格
    Mesh = dom.createElement('Mesh')
    Mesh.setAttribute('cmd', 'Import')
    root.appendChild(Mesh)
    ##
    MeshType = dom.createElement('MeshType')
    if meshtype == 1: text = dom.createTextNode('MID')
    if meshtype == 2: text = dom.createTextNode('3D')
    MeshType.appendChild(text)      # 网格类型: 壳体(MID), 实体(3D)
    Mesh.appendChild(MeshType)
    ##
    MeshUnit = dom.createElement('MeshUnit')
    text = dom.createTextNode('mm')
    MeshUnit.appendChild(text)      # 单位: mm
    Mesh.appendChild(MeshUnit)
    ##
    FileName = dom.createElement('FileName')
    text = dom.createTextNode(bdf_name)
    FileName.appendChild(text)      # 网格文件名称
    Mesh.appendChild(FileName)

    # 边界条件(注射口位置)
    BoundaryCondition = dom.createElement('BoundaryCondition')
    root.appendChild(BoundaryCondition)
    for i in range(len(inj_location)):
        ilocation = inj_location[i]
        ivector = inj_vector[i]
        ##
        InjLocation = dom.createElement('InjLocation')
        InjLocation.setAttribute('cmd', 'Create')
        BoundaryCondition.appendChild(InjLocation)
        ###
        if inj_type == 1: CoordinatesNormalized = dom.createElement('NodeID')
        if inj_type == 2: CoordinatesNormalized = dom.createElement('CoordinatesAbsolute')
        if inj_type == 3: CoordinatesNormalized = dom.createElement('CoordinatesNormalized')
        if inj_type == 1: text = dom.createTextNode('%d' % tuple(ilocation))
        if inj_type == 2: text = dom.createTextNode('%f %f %f' % tuple(ilocation))
        if inj_type == 3: text = dom.createTextNode('%f %f %f' % tuple(ilocation))
        CoordinatesNormalized.appendChild(text)      # 注射口位置  1:节点号  2:绝对坐标  3:归一化坐标
        InjLocation.appendChild(CoordinatesNormalized)
        ###
        TSetID = dom.createElement('TSetID')
        text = dom.createTextNode('40000')
        TSetID.appendChild(text)
        InjLocation.appendChild(TSetID)
        ###
        Vector = dom.createElement('Vector')
        text = dom.createTextNode('%f %f %f' % tuple(ivector))
        Vector.appendChild(text)
        InjLocation.appendChild(Vector)

    # 材料
    Material = dom.createElement('Material')
    Material.setAttribute('ID', str(mat_ID))
    Material.setAttribute('Shot', '1')
    root.appendChild(Material)

    # 求解器设置
    Property = dom.createElement('Property')
    root.appendChild(Property)
    ## 
    TSet = dom.createElement('TSet')
    Property.appendChild(TSet)
    ###
    ID = dom.createElement('ID')
    text = dom.createTextNode('10000')
    ID.appendChild(text)
    TSet.appendChild(ID)
    ###
    SubID = dom.createElement('SubID')
    text = dom.createTextNode('1')
    SubID.appendChild(text)
    TSet.appendChild(SubID)
    ###
    TCode = dom.createElement('TCode')
    TSet.appendChild(TCode)
    ####
    ID = dom.createElement('ID')
    text = dom.createTextNode('402')
    ID.appendChild(text)
    TCode.appendChild(ID)
    ####
    Value = dom.createElement('Value')
    text = dom.createTextNode('400')
    Value.appendChild(text)         # 最大熔体温度迭代次数
    TCode.appendChild(Value)
    ###
    TCode = dom.createElement('TCode')
    TSet.appendChild(TCode)
    ####
    ID = dom.createElement('ID')
    text = dom.createTextNode('52038')
    ID.appendChild(text)
    TCode.appendChild(ID)
    ####
    Value = dom.createElement('Value')
    text = dom.createTextNode('3')
    Value.appendChild(text)         # 自定义线程数
    TCode.appendChild(Value)
    ###
    TCode = dom.createElement('TCode')
    TSet.appendChild(TCode)
    ####
    ID = dom.createElement('ID')
    text = dom.createTextNode('52041')
    ID.appendChild(text)
    TCode.appendChild(ID)
    ####
    Value = dom.createElement('Value')
    text = dom.createTextNode(str(thread))
    Value.appendChild(text)         # 线程数
    TCode.appendChild(Value)

    f= open('study.xml', 'w')
    dom.writexml(f, addindent='    ', newl='\n')
    f.close()
    return

# /abaqus输出文件转化为moldflow计算文件
## ================================================== ##




## -------------------------------------------------- ##
# moldflow输出文件转化为abaqus输入文件

# 重命名文件
def RenameFile(oldname, newname):
    try:
        os.rename(oldname, newname)
    except WindowsError:
        os.remove(newname)
        os.rename(oldname, newname)

# 纤维主方向(单元局部坐标系)
def FiberPrincipalDirection(filename, NoE):

    f = open(filename,'r')
    rawdata = f.read()
    f.close()

    tree = ET.fromstring(rawdata)
    Dataset = tree.find('Dataset')
    Blocks = Dataset.find('Blocks')
    Block = Blocks[-1]
    Value = np.zeros([NoE, 6])
    EleID = list()

    # 数据储存至矩阵
    Data = Block[2]
    for i in range(len(Data)):
        ElementData = Data[i]
        ElementID = int(ElementData.attrib['ID'])  # 对应单元ID
        DeptValues = ElementData[0].text.split()
        for j in range(6):
            val = float(DeptValues[j])
            Value[ElementID-1,j] = val            # 数据存入 单元数*厚度方向数 的矩阵

    # 根据纤维取向张量计算主方向
    principal_direction_1st = np.zeros([NoE, 3])
    principal_direction_2nd = np.zeros([NoE, 3])
    for i in range(len(Value)):
        Vele = Value[i]
        mat = np.array([[Vele[0],     0.0,     0.0],
                        [Vele[3], Vele[1],     0.0],
                        [Vele[4], Vele[5], Vele[2]]])
        eigenvalue, featurevector = np.linalg.eigh(mat)

        principal_direction_1st[i] = featurevector[:,2]    # 第一主方向
        principal_direction_2nd[i] = featurevector[:,1]    # 第二主方向

    return principal_direction_1st, principal_direction_2nd

# 力学常数
def MechanicalConstant(filename, NoE):

    f = open(filename,'r')
    rawdata = f.read()
    f.close()

    tree = ET.fromstring(rawdata)
    Dataset = tree.find('Dataset')
    Blocks = Dataset.find('Blocks')
    Block = Blocks[-1]
    Value = np.zeros(NoE)

    # 数据储存至向量
    Data = Block[1]
    for i in range(len(Data)):
        ElementData = Data[i]
        ElementID = int(ElementData.attrib['ID'])
        DeptValues = ElementData[0].text.split()
        val = float(DeptValues[0])
        Value[ElementID-1] = val
    
    return Value

# 读取纤维取向结果和力学常数，写入inp附属文件
def InputIncludeFile(group, part, NoE):

    pd_1st, pd_2nd = FiberPrincipalDirection('group%d-fiber.xml' % group, NoE)
    E1  = MechanicalConstant('group%d-E1.xml' % group, NoE)
    E2  = MechanicalConstant('group%d-E2.xml' % group, NoE)
    G12 = MechanicalConstant('group%d-G12.xml' % group,NoE)
    Nu  = MechanicalConstant('group%d-Nu.xml' % group, NoE)

    # 单元Set
    filename = 'set.inp'
    f = open(filename,'w')
    for i in range(NoE):
        f.write('*Elset, elset=ELSET-%d\n %d,\n' % (i+1, i+1))
    f.close()

    # 材料属性
    filename = 'material.inp'
    f = open(filename,'w')
    for i in range(NoE):
        f.write('*Material, name=MATERIAL-%d\n' % (i+1))
        f.write('*Elastic, type=ENGINEERING CONSTANTS\n')
        f.write('%f, %f, %f, ' % (E1[i]*1e-6, E2[i]*1e-6, E2[i]*1e-6))
        f.write('%f, %f, %f, ' % (Nu[i], Nu[i], Nu[i]))
        f.write('%f, %f\n%f\n' % (G12[i]*1e-6, G12[i]*1e-6, G12[i]*1e-6))
    f.close()

    # 局部坐标系
    filename = 'orientation.inp'
    f = open(filename,'w')
    for i in range(NoE):
        f.write('*Orientation, name=ORI-%d\n' % (i+1))
        f.write('%f, %f, %f, %f, %f, %f\n' % tuple(list(pd_1st[i])+list(pd_2nd[i])))
        f.write('3, 0.\n')
    f.close()

    # 单元Section
    filename = 'section.inp'
    f = open(filename,'w')
    for section in part.Section:
        for el in part.ElementSet[section[1]]:
            f.write('*Shell Section, elset=ELSET-%d, material=MATERIAL-%d,' % (el, el))
            f.write(' orientation=ORI-%d\n%f, 5\n' % (el, section[3]))
    f.close()

    return

# 用于abaqus计算的inp文件
def NewInputFile(part, Others, filename):

    # 删除原先的文件
    for suffix in ['com', 'dat', 'inp', 'msg', 'odb', 'prt', 'sim', 'sta']:
        try: os.remove(filename[:-3]+suffix)
        except: pass
    
    f = open(filename,'w')

    # Part
    f.write('*Part, name=%s\n' %  part.Name)

    f.write('*Node\n')
    for i in range(len(part.Node)):
        f.write('%7d,%13f,%13f,%13f\n' % ((i+1), part.Node[i][0], part.Node[i][1], part.Node[i][2]))
    
    if 'S3' in part.ElementType.values():
        f.write('*Element, type=S3\n')
    for el in part.Element:
        element = part.Element[el]
        if part.ElementType[el] == 'S3':
            f.write('%7d,%7d,%7d,%7d\n' % (el, element[0], element[1], element[2]))
    
    f.write('*Include, input=orientation.inp\n')
    f.write('*Include, input=set.inp\n')
    f.write('*Include, input=section.inp\n')
    f.write('*End Part\n**\n')

    # Others
    for line in Others: 
        f.write(line)
        if '*End Assembly' in line:
            f.write('*Include, input=material.inp\n')
    
    f.close()
    return

# /moldflow输出文件转化为abaqus输入文件
## ================================================== ##



## -------------------------------------------------- ##
# 纤维取向张量和应力张量相似度计算

# 张量相似度计算
def mtx_similar(arr1:np.ndarray, arr2:np.ndarray, k:float, weighting:bool) ->float:

    # 求特征值和特征向量
    e1, v1 = np.linalg.eigh(arr1)
    e2, v2 = np.linalg.eigh(arr2)
    # 特征值归一化
    facter1 = np.linalg.norm(e1)
    facter2 = np.linalg.norm(e2)
    norm_e1 = abs(e1) / facter1
    norm_e2 = abs(e2) / facter2
    # 计算相似度
    similar = 0
    for i in range(3):
        for j in range(3):
            similar += abs(norm_e1[i]*norm_e2[j]) * (v1[:,i].dot(v2[:,j]))**2
    # 惩罚系数
    similar = np.power(similar, k)
    # 应力加权
    if weighting: similar *= facter1*facter2

    return similar

# 所有单元的相似度
def SimilarCalculate(NoE, fiberfile, stressfile, k:float):

    # 纤维取向张量
    f = open(fiberfile,'r')
    rawdata = f.read()
    f.close()
    tree = ET.fromstring(rawdata)
    Dataset = tree.find('Dataset')
    Blocks = Dataset.find('Blocks')
    Block = Blocks[-1]
    Fiber = np.zeros([NoE, 6])
    EleID = list()
    # 数据储存至矩阵
    Data = Block[2]
    for i in range(len(Data)):
        ElementData = Data[i]
        ElementID = int(ElementData.attrib['ID'])  # 对应单元ID
        DeptValues = ElementData[0].text.split()
        for j in range(6):
            val = float(DeptValues[j])
            Fiber[ElementID-1,j] = val            # 数据存入 单元数*厚度方向数 的矩阵

    # 应力张量
    Stress = np.loadtxt(stressfile)

    # 计算纤维取向张量和应力张量的相似度
    Elsimilar = np.zeros(NoE)
    Elsimilar_weighting = np.zeros(NoE)

    for i in range(NoE):

        Fele = Fiber[i]
        mat_fiber = np.array([[Fele[0], Fele[3], Fele[4]],
                              [Fele[3], Fele[1], Fele[5]],
                              [Fele[4], Fele[5], Fele[2]]])
        Sele = Stress[i]
        mat_stress = np.array([[Sele[0], Sele[3], Sele[4]],
                               [Sele[3], Sele[1], Sele[5]],
                               [Sele[4], Sele[5], Sele[2]]])
        # 计算单元相似度
        similar = mtx_similar(mat_fiber, mat_stress, k, False)
        similar_weighting = mtx_similar(mat_fiber, mat_stress, k, True)
        Elsimilar[i] = similar
        Elsimilar_weighting[i] = similar_weighting
    
    return Elsimilar, Elsimilar_weighting, Fiber

# 画相似度云图
def plot_contour(Coord, Value, figname):

    Coordx = Coord[:,0]
    Coordy = Coord[:,1]

    triang = tri.Triangulation(Coordx, Coordy)
    ##三角剖分
    fig = plt.figure(dpi=100, figsize=(8,6))
    ax = fig.subplots()
    # plt.tight_layout(pad = 0)
    plt.margins(0, 0)   # 使边框贴着图的边缘
    ax.set_aspect('equal')
    ax.set_xticks([])   #不显示x坐标轴
    ax.set_yticks([])   #不显示y坐标轴
    # im = ax.tricontourf(triang, Value, cmap='rainbow')
    Vmax = max(1, max(Value))
    im = ax.tripcolor(triang, Value, vmin=0.0, vmax=Vmax, cmap='Greys', shading='flat')
    # im = ax.tripcolor(triang, Value, cmap='Greys', shading='flat')
    cb = fig.colorbar(im)
    cb.ax.tick_params(labelsize=13)
    if figname: plt.savefig(figname)
    plt.show()


# /moldflow输出文件转化为abaqus输入文件
## ================================================== ##


if __name__ == "__main__":
    main()