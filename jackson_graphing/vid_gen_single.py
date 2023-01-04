import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
import os
from glob import glob
import argparse
color_list =["#3eb991","#e9a820","#e01563","#edb196","#6ecadc","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
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
    plt.savefig(f"../graphs/2D/{time_frame}")
    plt.close() 
if __name__ == '__main__':
    os.makedirs("../graphs/2D", exist_ok=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_name","-fn",type=str,help="the name of the output csv file")
    args = parser.parse_args()
    file_name = args.file_name
    propogate_ap_df = pd.DataFrame(pd.read_csv(f"../data_files/reduced2dV_output{file_name}.csv"))


    # removed the old graph png from the directory
    old_files = glob("../graphs/2D/*.png")
    for f in old_files:
        os.remove(f)

    for i in range(len(propogate_ap_df)): 
        graph_time_frame(propogate_ap_df,i)

    os.system(
                    "convert -delay 1 -loop 0 $(ls -1 ../graphs/2D/*.png | sort -V) -quality 95 ../vid/{}_{}.mp4".format("propagation",file_name))