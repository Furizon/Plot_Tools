import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
import matplotlib as mpl

def readFile(filePath):
    """把文件读取成列表形式

    Args:
        filePath (_type_): 文件路径

    Returns:
        list: 文件按行分隔的列表 
    """    

    fileHandler = open(filePath, "r")
    listOfLines = fileHandler.readlines()
    fileHandler.close()
    return listOfLines


def convertDSFfile(dsfPath):
    listOfLines = readFile(dsfPath)

    numOfLeads = int(listOfLines[1])
    typeOfLeads = listOfLines[3]
    
    # 读取原胞基矢
    unicellVectors = np.zeros((5, 3))
    for i in range(5, 10):
        unicellVectors[i - 5, :] = np.array(listOfLines[i].split(), dtype=float)

    numOfAtoms = int(listOfLines[11].split()[0])

    # 读取原子坐标
    atomsPosition = []
    for i in range(13, 13 + numOfAtoms):
        
        tup = (listOfLines[i].split()[0], np.array([eval(j) for j in listOfLines[i].split()[1:]]))
        atomsPosition.append(tup)

    # 读取k点数量
    kPoints = [int(i) for i in listOfLines[13 + numOfAtoms].split()]
    # 读取电荷密度
    electronDensity = []
    for i in tqdm(listOfLines[14 + numOfAtoms:]):
        for j in i.split():
            electronDensity.append(eval(j))
    electronDensity = np.array(electronDensity)
    electronDensity = np.array(electronDensity).reshape((kPoints[0], kPoints[1], kPoints[2], 2), order='F')
    
    return numOfLeads, typeOfLeads, unicellVectors, numOfAtoms, atomsPosition, kPoints, electronDensity

def showMultiLayerSpinPolarization(dsfPath):

    numOfLeads, typeOfLeads, unicellVectors, numOfAtoms, atomsPosition, kPoints, electronDensity = convertDSFfile(dsfPath)

    spinPolarization = np.squeeze(electronDensity[:, :, :, 0] - electronDensity[:, :, :, 1])
    print(np.mean(np.abs(spinPolarization)))
    # np.save("spinPolarization", spinPolarization)
    # spinPolarization = np.load("spinPolarization.npy")
    fig = plt.figure()
    ax = fig.add_subplot(121, projection='3d')
    # 绘制电荷密度分布
    X = np.linspace(0, unicellVectors[0, 0], kPoints[0])
    Y = np.linspace(0, unicellVectors[1, 1], kPoints[1])
    Y, X = np.meshgrid(Y, X)
    X = X + np.linspace(0, unicellVectors[1, 0], kPoints[0])
    Z = 36.18860614 * np.ones((kPoints[0], kPoints[1]))
    ax.scatter(X, Y, Z, c = spinPolarization[:, :, int((36.18860614 / unicellVectors[2, 2]) * kPoints[2])], zdir = 'z', s = 0.1, cmap='cool', depthshade=True)

    Z = 22.35136870 * np.ones((kPoints[0], kPoints[1]))
    ax.scatter(X, Y, Z, c = spinPolarization[:, :, int((22.35136870 / unicellVectors[2, 2]) * kPoints[2])], zdir = 'z', s = 0.1, cmap='cool', depthshade=True)

    Z = 29.26998763 * np.ones((kPoints[0], kPoints[1]))
    im = ax.scatter(X, Y, Z, c = spinPolarization[:, :, int((29.26998763 / unicellVectors[2, 2]) * kPoints[2])], zdir = 'z', s = 0.1, cmap='cool', depthshade=True)
    
    fig.colorbar(im)
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(X, Y, np.zeros((kPoints[0], kPoints[1])), c = spinPolarization[:, :, int((29.26998763 / unicellVectors[2, 2]) * kPoints[2])], cmap = 'cool', s = 1)
    ax2.view_init(90,-90)
    # 绘制原子
    elementTypes = ['Cr']
    elementColorsDict = {'Cr':'#000000', 'I':'#D65DEB'}
    for elementType in elementTypes:

        xs = [i[1][0] for i in atomsPosition if i[0] == elementType]
        ys = [i[1][1] for i in atomsPosition if i[0] == elementType]
        zs = [i[1][2] for i in atomsPosition if i[0] == elementType]
        ax.scatter(xs, ys, zs, c = elementColorsDict[elementType], marker = 'o', s = 3, alpha=1.0)

    
    plt.show()

def showLinearSpinPolarization(dsfPath):
    numOfLeads, typeOfLeads, unicellVectors, numOfAtoms, atomsPosition, kPoints, electronDensity = convertDSFfile(dsfPath)

    spinPolarization = np.squeeze(electronDensity[int(kPoints[0] / 2), int(kPoints[1] / 2), :, 0] - electronDensity[int(kPoints[0] / 2), int(kPoints[1] / 2), :, 1])
    fig = plt.figure()
    print(np.mean(np.abs(np.squeeze(electronDensity[:, :, :, 0] - electronDensity[:, :, :, 1]))))
    print(np.mean(np.squeeze(electronDensity[:, :, :, 0] - electronDensity[:, :, :, 1])))
    plt.plot(range(kPoints[2]), spinPolarization)
    plt.show()
    fig.savefig('svg/test.svg', transparent=True)

def plotTransmission(data, range):
    TransmissionCoefficients = np.fromstring(
        data.get("transm.transmissionCoefficients").replace("\n", ""),
        dtype=float,
        sep="          ",
    )
    EnergyPoints = np.fromstring(
        data.get("transm.energyPoints").replace("\n", ""), dtype=float, sep="          "
    )
    len = TransmissionCoefficients.shape[0] / 2
    TransmissionCoefficientsSpinUp = TransmissionCoefficients[:len]
    TransmissionCoefficientsSpinDown = TransmissionCoefficients[len:]
    plt.plot(EnergyPoints, TransmissionCoefficientsSpinUp[range], "r")
    plt.plot(EnergyPoints, TransmissionCoefficientsSpinDown[range], "b")
    plt.show()
if __name__ == '__main__':
    plotTransmission()
    # showLinearSpinPolarization("D:\Documents\Devicestudio\AgElectrodeCarbonChain\AgElectrodeCarbonChain\Ag_11C_Ag_h128_c200_Relaxed\Device\ElectronDensityTransmission\\2\TotalElectronDensity.dsf")
