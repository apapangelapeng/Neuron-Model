import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import os
import pathlibpt
from glob import glob
import glob

def graph_time_frame(data_frame,time_frame):
    row1 = data_frame.iloc[time_frame]
    x = np.arange(0,len(row1),step = 1)
    fig, ax = plt.subplots()
    ax.plot(x,row1)
    ax.set_xlabel("axon location")
    ax.set_ylabel("V")
    ax.set_title(f"time: {time_frame}")
    ax.set_ylim(-0.5,1.2)
    plt.plot()
    plt.tight_layout()
    plt.savefig(f"{current_path}/graphs/2D/{time_frame}")
    plt.close() 

    

current_path = str(pathlib.Path().resolve()).split("/")
print(current_path)
for index in range(len(current_path)):
    j = current_path[index]
    if j == "Neuron-Model":
        break
current_path ="/".join(current_path[:index+1])
print(current_path)
all_files = glob. glob(f"{current_path}/data_files/reduced2dV_output*.csv")
for f in all_files:
    propogate_ap_df = pd.DataFrame(pd.read_csv(f))
    name = ".".join(os.path.basename(f).split(".")[:-1])
    name ="_".join( name.split("_")[1:])
    print(name)
    #cleaning old files
    os.makedirs(f"{current_path}/graphs/2D", exist_ok=True)
    old_files =glob.glob(f"{current_path}/graphs/2D/*.png")
    for f in old_files:
        os.remove(f)
    print("removed old file")

    #graph new data
    for i in range(len(propogate_ap_df)): 
        graph_time_frame(propogate_ap_df,i)
    print("finished creating graphs, making video")

    os.system(f"convert -delay 1 -loop 0 $(ls -1 {current_path}/graphs/2D/*.png | sort -V) -quality 95 {current_path}/vid/{name}.mp4")