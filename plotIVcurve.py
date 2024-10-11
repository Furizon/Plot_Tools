import os
import re
import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.ticker import ScalarFormatter

def plotIVcurve(base_folder, prefix_pattern):
    # Initialize arrays to hold the data
    spinUpCurrent = []
    spinDownCurrent = []
    bias = []

    # Loop through each subfolder in the base folder
    for subfolder in os.listdir(base_folder):
        if os.path.isdir(os.path.join(base_folder, subfolder)) and re.match(prefix_pattern, subfolder):
            subfolder_path = os.path.join(base_folder, subfolder)
            

            
            # Construct the path to the CurrentVoltageCurves.xml file
            xml_file_path = os.path.join(subfolder_path, 'Device', 'CurrentVoltageCurves.xml')
            
            # Check if the XML file exists
            if os.path.isfile(xml_file_path):
                # Extract bias value from the folder name
                bias_value = subfolder.split('_')[-1]
                bias.append(float(bias_value))
                # Parse the XML file
                tree = ET.parse(xml_file_path)
                root = tree.getroot()
                
                # Find the ivc.I1_spinDecomposed data
                for elem in root.findall('.//nadata'):
                    name = elem.find('name').text
                    if name == 'ivc.I1_spinDecomposed':
                        values = elem[0].tail.strip().split()
                        spinUpCurrent.append(float(values[0]))
                        spinDownCurrent.append(float(values[1]))

    plt.rc('font',family='Times New Roman')
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.figure(figsize=(10, 6))
    ax = plt.axes(polar=False)
    # 设置 y 轴的显示格式
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((0, 0))
    ax.yaxis.get_offset_text().set_size(22)
    ax.yaxis.set_major_formatter(formatter)
    ax.plot(bias, spinUpCurrent, label='Spin Up', color="#cc3333", linestyle='-', linewidth=2, marker='+', markersize='10', markeredgewidth=1)
    ax.plot(bias, spinDownCurrent, label='Spin Down', color="#3333cc", linestyle='--', linewidth=2, marker='+', markersize='10', markeredgewidth=1)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5], fontsize=24)
    plt.yticks([0, 2.0e-6, 4.0e-6, 6.0e-6, 8.0e-6, 10.0e-6], fontsize=24)
    # plt.xlim(x_ticks[0] - 5, x_ticks[-1] + 5)
    plt.ylim(-5e-7, 1e-5)

    plt.xlabel("Bias(V)", fontsize=26)
    plt.ylabel('Current(A)', labelpad=26, rotation=90, fontsize=26)
    plt.tick_params(which='both', direction='in', length=4, width=1, colors='black')
    plt.legend(loc='best', fontsize=26, framealpha=0)
    plt.tight_layout()
    plt.show()

# Base folder where your subfolders are located
base_folder = 'D:\ProgramData\DeviceStudioProject\CarbonChainReplace\CarbonChainReplace'
prefix_pattern = r'^Ag_O_8C_O_Ag_h128_c200_Relaxed_([0-9]*\.?[0-9]+)$'

# base_folder = 'D:\ProgramData\DeviceStudioProject\AgElectrodeCarbonChain\AgElectrodeCarbonChain'
# prefix_pattern = r'^Ag_11C_Ag_h128_c200_Relaxed_([0-9]*\.?[0-9]+)$'

plotIVcurve(base_folder, prefix_pattern)