#!/usr/bin/env python
# coding: utf-8

# In[18]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import glob
color_list =["#3eb991","#6ecadc","#e9a820","#e01563","#edb196","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
test_output = pd.DataFrame(pd.read_csv("../data_files/test_output.csv"))

test_output = test_output.transpose()
test_output= test_output.reset_index(level=0)
print(test_output)
test_output = test_output.dropna()
column_name = list(test_output.iloc[0])
column_name = [ i.strip() for i in column_name]
test_output = test_output.drop(index=0)
test_output.columns = column_name
test_output["Tau_n"] = pd.to_numeric(test_output["Tau_n"])
print(test_output)


# # Data Processing

# # Graph 1
# (first 3)

# In[19]:


first_3 =list(test_output.columns )[:3]
fig, ax = plt.subplots()
x = np.arange(-40, 101, step=1)

for i in range(3):
    label = first_3[i]
    #ax.scatter(x,test_output[label],color = color_list[i])
    ax.plot(x,test_output[label],color = color_list[i])
    
colorlist = zip(first_3,color_list)
handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
t = ",".join(first_3)
plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
plt.plot()
ax.set_title(t)
plt.tight_layout()
plt.savefig(f"../graphs/tau/{t}.png")


# # Graph 2

# In[20]:


first_3 = list(test_output.columns)[3:3+3]
fig, ax = plt.subplots()
x = np.arange(-40, 101, step=1)
print(x.shape)

for i in range(3):
    label = first_3[i]
    #ax.scatter(x,test_output[label],color = color_list[i])
    ax.plot(x,test_output[label],color = color_list[i])
    
colorlist = zip(first_3,color_list)
handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
t = ",".join(first_3)
plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
plt.plot()
ax.set_title(t)
plt.tight_layout()
plt.savefig(f"../graphs/inf/{t}.png")


# # Graph 6 Voltages

# In[21]:


all_files = glob. glob(f"../data_files/testV_output_*.csv")
for f in all_files:
    voltage = pd.DataFrame(pd.read_csv(f))[:4]

    voltage = voltage.transpose().reset_index(level=0)
    name = str(os.path.basename(f).split(".")[0])
    name = name.split("_")[-1]
    name = f"(I_{name})"

    column =list(voltage.iloc[0])
    column = [ i.strip() for i in column]
    voltage = voltage.dropna()
    voltage = voltage.drop(index=0)
    voltage.columns = column
    issue_col = list(voltage.columns)[0]
    try:
        voltage[issue_col] = pd.to_numeric(voltage[issue_col])
    except:
        print(name)
    

    #print(voltage) # this shows the voltage 
    def graph_v_output(i, voltage):
        first_3 =  [list(voltage.columns)[i]]
        fig, ax = plt.subplots()
        size = len(voltage)
        x = np.arange(0,size, step=1)

        for i in range(len(first_3)):
            label = first_3[i]
            #ax.scatter(x,test_output[label],color = color_list[i])
            ax.plot(x,voltage[label],color = color_list[i])
            

        colorlist = zip(first_3,color_list)
        handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
        t = ",".join(first_3) + name
        ax.set_title(t)
        plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
        plt.plot()
        plt.tight_layout()
        plt.savefig(f"../graphs/voltage/{t}.png")

    graph_v_output(0, voltage)


    first_3 =  list(voltage.columns)[1:]
    fig, ax = plt.subplots()
    size = len(voltage)
    x = np.arange(0,size, step=1)

    for i in range(len(first_3)):
        label = first_3[i]
        ax.plot(x,voltage[label],color = color_list[i])
        

    colorlist = zip(first_3,color_list)
    handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
    t = ",".join(first_3) + name
    ax.set_title(t)
    plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
    plt.plot()
    plt.tight_layout()
    plt.savefig(f"../graphs/current/{t}.png")

