import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
class GraphStatic:
   
    def graph_df(data_frame:pd.DataFrame,begin_index=0,end_index=None):
        color_list =["#3eb991","#e9a820","#e01563","#edb196","#6ecadc","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
        if end_index:
            column_name =list(data_frame.columns )[begin_index:end_index]
        else:
            column_name =list(data_frame.columns )[begin_index:]
        fig, ax = plt.subplots()
        x = np.arange(0, len(data_frame)/1000, step=1/1000)

        for i in range(len(column_name)):
            label = column_name[i]
            ax.plot(x,data_frame[label],color = color_list[i])
            ax.set_xlabel('time (ms)')
            ax.set_ylabel('voltage (mV)')
            
        colorlist = zip(column_name,color_list)
        handles = [mpatches.Patch(color=colour, label=column_name) for label, colour in colorlist]
        t = ",".join(column_name)
        plt.legend(handles, column_name, ncol=1, bbox_to_anchor=(1, 1))
        plt.plot()
        ax.set_title(t)
        plt.tight_layout()

    def graph_df2(data_frame:pd.DataFrame,begin_index=0,end_index=None):
        color_list =["#3eb991","#e9a820","#e01563","#edb196","#6ecadc","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
        if end_index:
            column_name =list(data_frame.columns )[begin_index:end_index]
        else:
            column_name =list(data_frame.columns )[begin_index:]
        fig, ax = plt.subplots()
        x = np.arange(0, len(data_frame)/1000, step=1/1000)

        for i in range(len(column_name)):
            label = column_name[i]
            #ax.scatter(x,test_output[label],color = color_list[i])
            ax.plot(x,data_frame[label],color = color_list[i])
            ax.set_xlabel('time (ms)')
            ax.set_ylabel('current (uA)')
            
        colorlist = zip(column_name,color_list)
        handles = [mpatches.Patch(color=colour, label=column_name) for label, colour in colorlist]
        t = ",".join(column_name)
        plt.legend(handles, column_name, ncol=1, bbox_to_anchor=(1, 1))
        plt.plot()
        ax.set_title(t)
        plt.tight_layout()
