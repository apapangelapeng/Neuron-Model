#!/usr/bin/env python
# coding: utf-8

# In[161]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
color_list =["#3eb991","#6ecadc","#e9a820","#e01563","#edb196","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
test_output = pd.DataFrame(pd.read_csv("../data_files/test_output.csv"))

test_output = test_output.transpose()
test_output= test_output.reset_index(level=0)
test_output = test_output.dropna()
column_name = list(test_output.iloc[0])
column_name = [ i.strip() for i in column_name]
test_output = test_output.drop(index=0)
test_output.columns = column_name
test_output["N"] = pd.to_numeric(test_output["N"])



# # Data Processing

# # Graph 1
# (first 3)

# In[162]:


first_3 =list(test_output.columns )[:3]
fig, ax = plt.subplots()
x = np.arange(-40, 101, step=1)
for i in range(len(first_3)):
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
plt.savefig(f"../graphs/{t}.png")


# # Graph 2

# In[163]:


first_3 = list(test_output.columns)[3:3+3]
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
plt.savefig(f"../graphs/{t}.png")
plt.show()

# # Graph3

# In[164]:


first_3 = list(test_output.columns)[6:]
fig, ax = plt.subplots()
x = np.arange(-40, 101, step=1)


for i in range(2):
    label = first_3[i]
    #ax.scatter(x,test_output[label],color = color_list[i])
    ax.plot(x,test_output[label],color = color_list[i])
    
colorlist = zip(first_3,color_list)
handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
t = ",".join(first_3)
ax.set_title(t)
plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
plt.plot()
plt.tight_layout()
plt.savefig(f"../graphs/{t}.png")
plt.show()
