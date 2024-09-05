import xml.etree.ElementTree as ET
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import os
from tqdm import tqdm
import plotly.graph_objects as go
from matplotlib.ticker import ScalarFormatter
###############玻尔单位制换算##############
# 1 a.u.（长度）= a0 = 0.5291772083×10-10 m = 0.5291772083 Å
BOHR_LENGTH = 0.5291772083
HATREE = 27.2114
def readDSFfile(FilePath):
    Data = {}
    try:
        FileHandler = open(FilePath, "r")
        ListOfLines = FileHandler.readlines()
    except Exception as e:
        print("Error occured while reading files %s" % (e))
    finally:
        FileHandler.close()

    # 电极信息
    NumOfLeads = int(ListOfLines[1])
    TypeOfLeads = ListOfLines[3]

    # 读取原胞基矢
    UnicellVectors = np.zeros((5, 3))
    for i in range(5, 9):
        UnicellVectors[i - 5, :] = np.array(ListOfLines[i].split(), dtype=float)

    # 原子数量
    NumOfAtoms = int(ListOfLines[11].split()[0])

    # 读取原子坐标
    AtomPositions = []
    for i in range(13, 13 + NumOfAtoms):
        tup = (
            ListOfLines[i].split()[0],
            np.array([eval(j) for j in ListOfLines[i].split()[1:]]),
        )
        AtomPositions.append(tup)

    # 读取空间采样点数量 K点
    KPoints = [int(i) for i in ListOfLines[13 + NumOfAtoms].split()]
    # 读取电荷密度
    ElectronDensity = []
    for i in tqdm(ListOfLines[14 + NumOfAtoms :]):
        for j in i.split():
            ElectronDensity.append(eval(j))
    ElectronDensity = np.array(ElectronDensity)
    ElectronDensity = np.array(ElectronDensity).reshape(
        (KPoints[0], KPoints[1], KPoints[2], 2), order="F"
    )

    Data["NumOfLeads"] = NumOfLeads
    Data["TypeOfLeads"] = TypeOfLeads
    Data["UnicellVectors"] = UnicellVectors
    Data["NumOfAtoms"] = NumOfAtoms
    Data["AtomPositions"] = AtomPositions
    Data["KPoints"] = KPoints
    Data["ElectronDensity"] = ElectronDensity
    return Data

def parseElement(root):
    if root.find("nadata") == None:
        if len(list(root)) == 0:
            return {root.text: ""}
        if len(list(root)) == 1:
            return {root[0].text: root[0].tail}
        elif len(list(root)) == 2:
            return {root[0].text: root[1].text}

    else:
        temp = {}
        for child in list(root):
            temp.update(parseElement(child))
        return temp

def plotTransmission(Path):
    xml_data = open(Path, "r").read()  # Read file
    Data = parseElement(ET.XML(xml_data))
    TransmissionCoefficients = np.fromstring(
        Data.get("transm.transmissionCoefficients").replace("\n", ""),
        dtype=float,
        sep="          ",
    )
    EnergyPoints = np.fromstring(
        Data.get("transm.energyPoints").replace("\n", ""), dtype=float, sep="          "
    ) * HATREE
    plt.figure(figsize=(8, 6))
    ax = plt.axes(polar=False)
    
    len = int(TransmissionCoefficients.shape[0] / 2)
    TransmissionCoefficientsSpinUp = TransmissionCoefficients[:len]
    TransmissionCoefficientsSpinDown = TransmissionCoefficients[len:]
    
    ax.plot(EnergyPoints, TransmissionCoefficientsSpinUp, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
    ax.plot(EnergyPoints, TransmissionCoefficientsSpinDown, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks([-1.0, -0.7, -0.4, 0, 0.4, 0.7, 1.0], fontsize=22)
    plt.yticks(np.linspace(0, 3, 7), fontsize=22)
    plt.xlim(-1, 1)
    # plt.ylim(0, 4)
    plt.xlabel(r'$E-E_{F}(eV)$', fontsize=40)
    ax.set_ylabel('Transmission', labelpad=15, fontsize=40)
    ax.spines['top'].set_color('black')    # 设置顶部边框线颜色为灰色
    ax.spines['right'].set_color('black')  # 设置右侧边框线颜色为灰色
    ax.spines['bottom'].set_color('black') # 设置底部边框线颜色为黑色
    ax.spines['left'].set_color('black')  # 设置左侧边框线颜色为黑色
    ax.patch.set_alpha(0.0)   # 轴区域的背景
    ax.tick_params(which='both', direction='in', length=4, width=1, colors='black')
    plt.legend(loc='upper right', fontsize=30, framealpha=0)
    plt.tight_layout()
    plt.show()

def calSpinSpilting(data):
    TransmissionCoefficients = np.fromstring(
        data.get("transm.transmissionCoefficients").replace("\n", ""),
        dtype=float,
        sep="          ",
    )
    len = int(TransmissionCoefficients.shape[0] / 2)
    TransmissionCoefficientsSpinUp = TransmissionCoefficients[:len]
    TransmissionCoefficientsSpinDown = TransmissionCoefficients[len:]
    EnergyPoints = np.fromstring(
        data.get("transm.energyPoints").replace("\n", ""), dtype=float, sep="          "
    )
    SpinSpilting = (TransmissionCoefficientsSpinUp[201] - TransmissionCoefficientsSpinDown[201]) / (TransmissionCoefficientsSpinUp[201] + TransmissionCoefficientsSpinDown[201])
    Transmission = TransmissionCoefficientsSpinUp[201]
    return SpinSpilting, Transmission

def readDFTBLocalCurrentFiles(FilePath):
    
    # 存储原子坐标
    AtomPositions = []
    # 读入supercell文件
    try:
        FileHandler = open(FilePath + 'supercell.xyz', "r")
        ListOfLines = FileHandler.readlines()
    except Exception as e:
        print("Error occured while reading files %s" % (e))
    finally:
        FileHandler.close()
    
    AtomsNum = eval(ListOfLines[0])
    for i in range(2, AtomsNum + 2):
        tup = (ListOfLines[i].split()[0], np.array([eval(j) for j in ListOfLines[i].split()[1:]]))
        AtomPositions.append(tup)
        
        
    # 读入lcurrents文件
    try:
        FileHandler = open(FilePath + 'lcurrents_u.dat', "r")
        ListOfLines = FileHandler.readlines()
    except Exception as e:
        print("Error occured while reading files %s" % (e))
    finally:
        FileHandler.close()
    
    
    # 建立一个矩阵，以存储原子间电流
    NeighboursCurrent = np.zeros((AtomsNum, AtomsNum))
    # 进行遍历
    for Line in ListOfLines:
        Line = Line.split()
        AtomNo = eval(Line[0])
        for i in range(6, len(Line), 2):
            NeighboursCurrent[AtomNo - 1, int(Line[i - 1]) - 1] = eval(Line[i])
    
    return AtomPositions, np.matrix(NeighboursCurrent)


def plotElectronDistribution(DSFdata):

    x_ticks = []
    NonAgAtomNum = 0
    for atom in DSFdata['AtomPositions']:
        if atom[0] != 'Ag':
            x_ticks.append((int)(atom[1][2] / DSFdata['UnicellVectors'][2][2] * DSFdata['KPoints'][2]))
            NonAgAtomNum += 1
    x_ticks = np.sort(x_ticks)
    x_labels = np.arange(1, NonAgAtomNum + 1, dtype=int)
    SpinUpDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, :, 0], axis = (0, 1))
    SpinDownDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, :, 1], axis = (0, 1))
    spinUpRatio = SpinUpDensitySumOverXY / (
        np.sum(SpinUpDensitySumOverXY[x_ticks[0] - 5:x_ticks[-1] + 5] + SpinDownDensitySumOverXY[x_ticks[0] - 5:x_ticks[-1] + 5])
    )
    spinDownRatio = SpinDownDensitySumOverXY / (
        np.sum(SpinUpDensitySumOverXY[x_ticks[0] - 5:x_ticks[-1] + 5] + SpinDownDensitySumOverXY[x_ticks[0] - 5:x_ticks[-1] + 5])
    )
    plt.figure(figsize=(8, 6))
    ax = plt.axes(polar=False)
    # 设置 y 轴的显示格式
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((0, 0))
    ax.yaxis.get_offset_text().set_size(22)
    ax.yaxis.set_major_formatter(formatter)
    ax.plot(range(spinUpRatio.shape[0]), spinUpRatio, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
    ax.plot(range(spinDownRatio.shape[0]), spinDownRatio, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
    ax.patch.set_alpha(0.0)   # 轴区域的背景
    
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(x_ticks, x_labels, fontsize=22)
    plt.yticks(fontsize=22)
    
    plt.xlim(x_ticks[0] - 5, x_ticks[-1] + 5)
    plt.ylim(0.002, 0.015)
    plt.xlabel("Sites", fontsize=40)
    ax.set_ylabel(r'$\rho$', labelpad=20, rotation=0, fontsize=40)
    ax.spines['top'].set_color('black')    # 设置顶部边框线颜色为灰色
    ax.spines['right'].set_color('black')  # 设置右侧边框线颜色为灰色
    ax.spines['bottom'].set_color('black') # 设置底部边框线颜色为黑色
    ax.spines['left'].set_color('black')  # 设置左侧边框线颜色为黑色
    ax.tick_params(which='both', direction='in', length=4, width=1, colors='black')
    plt.legend(loc='best', fontsize=30, framealpha=0)
    plt.tight_layout()
    plt.show()
    SpinPolarization = np.average(np.abs(spinUpRatio - spinDownRatio))
    print(SpinPolarization)
    
    
def cal():
    RootPath = (
        "D:\ProgramData\DeviceStudioProject\CarbonChainReplace\CarbonChainReplace"
    )
    Files = os.listdir(RootPath)
    for File in Files:
        if os.path.isdir(RootPath + "/" + File):
            Path = RootPath + "/" + File + "/Device/Transmission.xml"
            if os.path.isfile(Path):
                xml_data = open(Path, "r").read()  # Read file
                Data = parseElement(ET.XML(xml_data))
                print(File)
                print("{:.5f}".format(calSpinSpilting(Data, [0, 400])))
            Path = RootPath + "/" + File + "/Device/TotalElectronDensity.dsf"
            if os.path.isfile(Path):
               plotElectronDistribution(Path)

def plotPolarization(DSFdata):

    x_ticks = []
    NonAgAtomNum = 0
    for atom in DSFdata['AtomPositions']:
        if atom[0] != 'Ag':
            x_ticks.append((int)(atom[1][2] / DSFdata['UnicellVectors'][2][2] * DSFdata['KPoints'][2]))
            NonAgAtomNum += 1
    x_ticks = np.sort(x_ticks)
    x_labels = np.arange(1, NonAgAtomNum + 1, dtype=int)
    SpinUpDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, :, 0], axis = (0, 1))
    SpinDownDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, :, 1], axis = (0, 1))
    spinUpRatio = SpinUpDensitySumOverXY / (
        SpinUpDensitySumOverXY + SpinDownDensitySumOverXY
    )
    spinDownRatio = SpinDownDensitySumOverXY / (
        SpinUpDensitySumOverXY + SpinDownDensitySumOverXY
    )
    plt.figure(figsize=(8, 6))
    ax = plt.axes(polar=False)
    ax.plot(range(spinUpRatio.shape[0]), spinUpRatio, color="#ff7f0e", linestyle='-', linewidth=2, label="Spin Up")
    ax.plot(range(spinDownRatio.shape[0]), spinDownRatio, color="#1f77b4", linestyle='--', linewidth=2, label="Spin Down")
    
    
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(x_ticks, x_labels, fontsize=24)
    plt.yticks([0.3, 0.4, 0.5, 0.6, 0.7], fontsize=24)
    
    plt.xlim(x_ticks[0] - 5, x_ticks[-1] + 5)
    plt.ylim(0.4, 0.6)
    plt.xlabel("Sites", fontsize=26)
    ax.set_ylabel(r'$n_{\mathrm{i}}$', labelpad=20, rotation=0, fontsize=26)
    ax.spines['top'].set_color('black')    # 设置顶部边框线颜色为灰色
    ax.spines['right'].set_color('black')  # 设置右侧边框线颜色为灰色
    ax.spines['bottom'].set_color('black') # 设置底部边框线颜色为黑色
    ax.spines['left'].set_color('black')  # 设置左侧边框线颜色为黑色
    ax.tick_params(which='both', direction='in', length=4, width=1, colors='black')
    plt.legend(loc='best', fontsize=26, framealpha=0)
    plt.tight_layout()
    plt.show()
    SpinPolarization = np.average(np.abs(spinUpRatio - spinDownRatio))
    print(SpinPolarization)
    
def plotSpinSpiltingBias(RootPath):
    Files = os.listdir(RootPath)
    Bias = []
    Polarizations = []
    Transmissions = []
    for File in Files:
        if os.path.isdir(RootPath + "/" + File):
            Path = RootPath + "/" + File + "/Transmission.xml"
            if os.path.isfile(Path):
                xml_data = open(Path, "r").read()  # Read file
                Data = parseElement(ET.XML(xml_data))
                Bias.append(float(File))
                Polarization, Transmission = calSpinSpilting(Data)
                Polarizations.append(Polarization)
                Transmissions.append(Transmission)
    return Bias, Polarizations, Transmissions
def plotPolarizationBias(RootPath):
    Files = os.listdir(RootPath)
    Bias = []
    SpinSpilting = []
    for File in Files:
        if os.path.isdir(RootPath + "/" + File):
            Path = RootPath + "/" + File + "/TotalElectronDensity.dsf"
            if os.path.isfile(Path):
                DSFdata = readDSFfile(Path)
                NonAgAtomNum = 0
                x_ticks = []
                for atom in DSFdata['AtomPositions']:
                    if atom[0] != 'Ag':
                        x_ticks.append((int)(atom[1][2] / DSFdata['UnicellVectors'][2][2] * DSFdata['KPoints'][2]))
                        NonAgAtomNum += 1
                x_ticks = np.sort(x_ticks)
                SpinUpDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, x_ticks[0]:x_ticks[-1], 0], axis = (0, 1))
                SpinDownDensitySumOverXY = np.sum(DSFdata["ElectronDensity"][:, :, x_ticks[0]:x_ticks[-1], 1], axis = (0, 1))
                spinUpRatio = SpinUpDensitySumOverXY / (
                    SpinUpDensitySumOverXY + SpinDownDensitySumOverXY
                )
                spinDownRatio = SpinDownDensitySumOverXY / (
                    SpinUpDensitySumOverXY + SpinDownDensitySumOverXY
                )
                Bias.append(float(File))
                SpinSpilting.append(np.average(np.abs(spinUpRatio - spinDownRatio)))
    return Bias, SpinSpilting
# 判断两个原子是否成键
def is_bonded(atom1, atom2, bond_lengths):
    # 查找键长
    key = f"{atom1[0]}-{atom2[0]}"
    reverse_key = f"{atom2[0]}-{atom1[0]}"
    if key in bond_lengths:
        max_length = bond_lengths[key]
    elif reverse_key in bond_lengths:
        max_length = bond_lengths[reverse_key]
    else:
        return False
    # 判断距离
    distance = np.sqrt((atom1[1][0] - atom2[1][0]) ** 2 + (atom1[1][1] - atom2[1][1]) ** 2 + (atom1[1][2] - atom2[1][2]) ** 2)
    return distance <= max_length

# 生成球体的函数
def create_sphere(center, radius=0.5, resolution=10):
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    return x, y, z

# 生成圆柱体的函数
def create_cylinder(start, end, radius=0.2, resolution=10):
    v = np.array(end) - np.array(start)
    mag = np.linalg.norm(v)
    v = v / mag
    not_v = np.array([1, 0, 0]) if (v == np.array([0, 1, 0])).all() else np.array([0, 1, 0])
    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)
    t = np.linspace(0, mag, resolution)
    theta = np.linspace(0, 2 * np.pi, resolution)
    t, theta = np.meshgrid(t, theta)
    X, Y, Z = [
        start[i] + v[i] * t + radius * np.sin(theta) * n1[i] + radius * np.cos(theta) * n2[i]
        for i in [0, 1, 2]
    ]
    return X, Y, Z
def plotIsosurfaceSpinPolarization(DSFdata, iso_value):

    densityUp = DSFdata["ElectronDensity"][:, :, :, 0] - DSFdata["ElectronDensity"][:, :, :, 1]
    densityDown = DSFdata["ElectronDensity"][:, :, :, 1] - DSFdata["ElectronDensity"][:, :, :, 0]
    # 计算步长
    steps = [DSFdata["UnicellVectors"][i][i] / DSFdata["KPoints"][i] for i in range(3)]
    X, Y, Z = np.mgrid[0:DSFdata["UnicellVectors"][0][0]:steps[0], 
                       0:DSFdata["UnicellVectors"][1][1]:steps[1],
                       0:DSFdata["UnicellVectors"][2][2]:steps[2]]
    # 创建图形
    fig = go.Figure()
    
    
    # 绘制原子
    atoms = DSFdata['AtomPositions']
    for atom in atoms:
        x, y, z = create_sphere((atom[1][0], atom[1][1], atom[1][2]))
        fig.add_trace(go.Surface(x=x, y=y, z=z, colorbar=None, showscale=False))

    # 键长表
    bond_lengths = {
        'C-C': 1.54,
        'C-Ag': 2.04,
        'Ag-Ag': 2.88
    }
    for i, atom1 in enumerate(atoms):
        for j, atom2 in enumerate(atoms):
            if i < j and is_bonded(atom1, atom2, bond_lengths):
                x, y, z = create_cylinder((atom1[1][0], atom1[1][1], atom1[1][2]), (atom2[1][0], atom2[1][1], atom2[1][2]))
                fig.add_trace(go.Surface(x=x, y=y, z=z, colorbar=None, showscale=False))


    # 添加第一组数据的等能面（红色）
    fig.add_trace(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=densityUp.flatten(),
        isomin=iso_value,
        isomax=iso_value,
        surface_count=1,
        colorscale='Reds',
        opacity=1,
        caps=dict(x_show=False, y_show=False)
    ))

    # 添加第二组数据的等能面（蓝色）
    fig.add_trace(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=densityDown.flatten(),
        isomin=iso_value,
        isomax=iso_value,
        surface_count=1,
        colorscale='Blues',
        opacity=1,
        caps=dict(x_show=False, y_show=False)
    ))

    # 设置图表布局
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title='',
                showgrid=False,
                zeroline=False,
                visible=False,
                backgroundcolor='rgba(0,0,0,0)'
            ),
            yaxis=dict(
                title='',
                showgrid=False,
                zeroline=False,
                visible=False,
                backgroundcolor='rgba(0,0,0,0)'
            ),
            zaxis=dict(
                title='',
                showgrid=False,
                zeroline=False,
                visible=False,
                backgroundcolor='rgba(0,0,0,0)'
            ),

            camera=dict(
                eye=dict(x=0., y=2.5, z=0.),
                projection=dict(type='orthographic')
            )
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        width=3600,  # 宽度（像素）
        height=2400,  # 高度（像素）
        title='3D Isosurfaces with Atoms and Bonds (Orthographic Projection)',
        showlegend=False
    )

    # 显示图表
    fig.show()


if __name__ == "__main__":
    plt.rc('font',family='Times New Roman')
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    
    # Path = 'D:\ProgramData\DeviceStudioProject\AgElectrodeCarbonChain\AgElectrodeCarbonChain\Ag_11C_Ag_h128_c200_Relaxed\Device\\'
    Path = 'D:\ProgramData\DeviceStudioProject\AgElectrodeCarbonChain\AgElectrodeCarbonChain\Ag_11C_Ag_h128_c200_Relaxed\Device\\'
    # Path = 'D:\ProgramData\DeviceStudioProject\AgElectrodeCarbonChain\AgElectrodeCarbonChain\Ag_11C_Ag_h128_c100_Relaxed\Device\\'
    DSFdata = readDSFfile(Path + "EffectivePotential.dsf")
    # plotIsosurfaceSpinPolarization(DSFdata, iso_value = 0.01)
    plotElectronDistribution(DSFdata)
    # plotTransmission(Path + 'Transmission.xml')
    
    # Bias, Polarizations, Transmissions = plotSpinSpiltingBias('D:\ProgramData\DeviceStudioProject\data')
    # plt.figure(figsize=(8, 6))
    # ax = plt.axes(polar=False)
    # ax.plot(Bias, Polarizations, color='black', marker='+'
    #      ,markeredgecolor='black',markersize='10',markeredgewidth=1)
    # plt.grid(True, linestyle='--', alpha=0.7)

    # plt.xlabel("Bias(V)", fontsize=26)
    # ax.set_ylabel('Spin Filtering Ratio', labelpad=15, rotation=90, fontsize=26)
    # ax.spines['top'].set_color('black')    # 设置顶部边框线颜色为灰色
    # ax.spines['right'].set_color('black')  # 设置右侧边框线颜色为灰色
    # ax.spines['bottom'].set_color('black') # 设置底部边框线颜色为黑色
    # ax.spines['left'].set_color('black')  # 设置左侧边框线颜色为黑色
    # ax.tick_params( length=4, width=1, colors='black', labelsize=22)
    # ax.set_yticks(np.linspace(-1, 1, 5))
    # ax.set_ylim(-1.1, 1.2)
    # ax2 = ax.twinx()
    # ax2.set_ylabel('Transmission', labelpad=15, rotation=90, fontsize=26, color='#c94737', )
    # # 设置 y 轴的显示格式
    # formatter = ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True) 
    # formatter.set_powerlimits((0, 0))
    # ax2.yaxis.get_offset_text().set_size(22)
    # ax2.yaxis.set_major_formatter(formatter)
    # ax2.plot(Bias, Transmissions, color='#c94737', marker='+'
    #      ,markeredgecolor='#c94737',markersize='10',markeredgewidth=1)
    # ax2.tick_params(length=4, width=1, colors='#c94737', labelsize=22)
    # ax2.set_yticks(np.linspace(0, 2, 5))
    # ax2.set_ylim(-0.1, 2.2)
    
    # plt.tight_layout()
    # plt.show()
