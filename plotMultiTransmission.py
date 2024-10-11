import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
import os
import matplotlib
plt.rc('font',family='Times New Roman')
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['mathtext.default'] = 'regular'
###############玻尔单位制换算##############
# 1 a.u.（长度）= a0 = 0.5291772083×10-10 m = 0.5291772083 Å
BOHR_LENGTH = 0.5291772083
HATREE = 27.2114
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
def getTransmission(Path):
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
    
    len = int(TransmissionCoefficients.shape[0] / 2)
    TransmissionCoefficientsSpinUp = TransmissionCoefficients[:len]
    TransmissionCoefficientsSpinDown = TransmissionCoefficients[len:]
    
    return EnergyPoints, TransmissionCoefficientsSpinUp, TransmissionCoefficientsSpinDown

# 创建一些示例数据
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.tan(x)
y4 = np.exp(-x)

# 创建4个纵向排列的子图，并共享一个x轴
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(8, 8))

# 绘制每个子图
EnergyPoints, TransmissionCoefficientsSpinUp, TransmissionCoefficientsSpinDown = getTransmission(
    "D:\ProgramData\DeviceStudioProject\data\\0\Transmission.xml")
axs[0].plot(EnergyPoints, TransmissionCoefficientsSpinUp, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
axs[0].plot(EnergyPoints, TransmissionCoefficientsSpinDown, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
axs[0].legend(loc='upper right', fontsize=15, framealpha=0)
axs[0].grid(True, linestyle='--', alpha=0.7)
axs[0].set_ylim(-0.05, 2.1)
axs[0].set_xlim(-1, 1)
axs[0].text(1.02, 0.5, "Bias = 0V", transform=axs[0].transAxes, fontsize=18, verticalalignment='center')
axs[0].tick_params(axis='y', labelsize=15)

EnergyPoints, TransmissionCoefficientsSpinUp, TransmissionCoefficientsSpinDown = getTransmission(
    "D:\ProgramData\DeviceStudioProject\data\\0.2\Transmission.xml")
axs[1].plot(EnergyPoints, TransmissionCoefficientsSpinUp, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
axs[1].plot(EnergyPoints, TransmissionCoefficientsSpinDown, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
axs[1].grid(True, linestyle='--', alpha=0.7)
axs[1].set_ylim(-0.05, 2.1)
axs[1].set_xlim(-1, 1)
axs[1].text(1.02, 0.5, "           0.2V", transform=axs[1].transAxes, fontsize=18, verticalalignment='center')
axs[1].tick_params(axis='y', labelsize=15)

EnergyPoints, TransmissionCoefficientsSpinUp, TransmissionCoefficientsSpinDown = getTransmission(
    "D:\ProgramData\DeviceStudioProject\data\\0.25\Transmission.xml")
axs[2].plot(EnergyPoints, TransmissionCoefficientsSpinUp, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
axs[2].plot(EnergyPoints, TransmissionCoefficientsSpinDown, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
axs[2].grid(True, linestyle='--', alpha=0.7)
axs[2].set_ylim(-0.05, 2.1)
axs[2].set_xlim(-1, 1)
axs[2].text(1.02, 0.5, "           0.25V", transform=axs[2].transAxes, fontsize=18, verticalalignment='center')
axs[2].tick_params(axis='y', labelsize=15)

EnergyPoints, TransmissionCoefficientsSpinUp, TransmissionCoefficientsSpinDown = getTransmission(
    "D:\ProgramData\DeviceStudioProject\data\\0.35\Transmission.xml")
axs[3].plot(EnergyPoints, TransmissionCoefficientsSpinUp, color="#cc3333", linestyle='-', linewidth=2, label="Spin Up")
axs[3].plot(EnergyPoints, TransmissionCoefficientsSpinDown, color="#3333cc", linestyle='--', linewidth=2, label="Spin Down")
axs[3].grid(True, linestyle='--', alpha=0.7)
axs[3].set_ylim(-0.05, 2.1)
axs[3].set_xlim(-1, 1)
axs[3].text(1.02, 0.5, "           0.35V", transform=axs[3].transAxes, fontsize=18, verticalalignment='center')
axs[3].tick_params(axis='y', labelsize=15)
axs[3].tick_params(axis='x', labelsize=15)
# 设置x轴的标签
axs[3].set_xlabel(r'$E-E_{F}(eV)$', fontsize=18)
fig.text(0.04, 0.5, 'Transmission', va='center', rotation='vertical', fontsize=18)
# 使用subplots_adjust调整间距
plt.subplots_adjust(left=0.1, right=0.82, top=0.95, bottom=0.1, hspace=0.1)

# 显示图像
plt.show()
