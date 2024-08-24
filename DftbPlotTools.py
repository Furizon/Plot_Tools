import numpy as np
from matplotlib import pyplot as plt
class DftbOutput:
    
    # 读取电荷文件
    def readChargeFile(filePath):
        fileHandler = open(filePath, "r")
        listOfLines = fileHandler.readlines()
        fileHandler.close()
        chargeFileConverted = {}
        chargeFileConverted['version'] = int(listOfLines[0])

        informationLine = listOfLines[1]
        informationLine = informationLine.split()
        chargeFileConverted['flags'] = informationLine[0:4]
        chargeFileConverted['numOfSites'] = int(informationLine[4])
        chargeFileConverted['nSpin'] = int(informationLine[5])

        # 校验和
        chargeFileConverted['chargeCheckSum'] = eval(informationLine[6])
        chargeFileConverted['magnetisationCheckSum'] = eval(informationLine[7])

        # 读取电荷矩阵
        charge = []
        for line in listOfLines[2 : 2 + chargeFileConverted['numOfSites']]:
            for j in line.split():
                charge.append(eval(j))
        charge = np.array(charge)
        charge = charge.reshape(chargeFileConverted['numOfSites'], 4)
        chargeFileConverted['charge'] = charge


        # 读取磁矩矩阵
        magnetisation = []
        for line in listOfLines[2 + chargeFileConverted['numOfSites'] : 2 + 2 * chargeFileConverted['numOfSites']]:
            for j in line.split():
                magnetisation.append(eval(j))
        magnetisation = np.array(magnetisation)
        magnetisation = magnetisation.reshape(chargeFileConverted['numOfSites'], 4)
        chargeFileConverted['magnetisation'] = magnetisation
        if np.isclose(np.sum(chargeFileConverted['charge']), chargeFileConverted['chargeCheckSum']):
            if np.isclose(np.sum(chargeFileConverted['magnetisation']), chargeFileConverted['magnetisationCheckSum']):
                return chargeFileConverted
        
        print("CheckSum does not coherent with the data")
        return chargeFileConverted
    
    def calTotalSpinPolarization(chargeFileDict):
        return np.mean(np.abs(chargeFileConverted['magnetisation'][:, 0] - chargeFileConverted['magnetisation'][:, 1]))




if __name__ == '__main__':
    X = []
    Y = []
    for i in range(40):
        X.append(i)
        chargeFileConverted = DftbOutput.readChargeFile('C:\\Users\Admin\Downloads\charges\\' + str(i))
        Y.append(DftbOutput.calTotalSpinPolarization(chargeFileConverted))

    plt.plot(X, Y)
    plt.xlabel('Bias Voltage(100mV)')
    plt.ylabel('SpinPolarisation')
    plt.show()