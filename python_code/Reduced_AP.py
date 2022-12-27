#!/usr/bin/env python
# coding: utf-8

# In[52]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
color_list =["#3eb991","#6ecadc","#e9a820","#e01563","#edb196","#1f94ac","#ae9a6a","#ccb8a6","#343a44"]
test_output = pd.DataFrame(pd.read_csv("../data_files/reduced.csv"))

test_output = test_output.transpose()
test_output= test_output.reset_index(level=0)
print(test_output)
test_output = test_output.dropna()
column_name = list(test_output.iloc[0])
column_name = [ i.strip() for i in column_name]
test_output = test_output.drop(index=0)
test_output.columns = column_name
test_output["N"] = pd.to_numeric(test_output["N"])
print(test_output)


# # Data Processing

# # Graph 1
# (first 3)

# In[53]:


first_3 =list(test_output.columns )[:3]
fig, ax = plt.subplots()
x = np.arange(-40, 98, step=1)
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

# In[54]:


first_3 = list(test_output.columns)[3:3+3]
fig, ax = plt.subplots()
x = np.arange(-40, 98, step=1)
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
plt.savefig(f"../graphs/{t}.png")


# # Graph 3

# In[55]:


first_3 = list(test_output.columns)[6:6+3]
fig, ax = plt.subplots()
x = np.arange(-40, 98, step=1)
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
plt.savefig(f"../graphs/{t}.png")


# # Graph4 Proportions

# In[56]:


first_3 = list(test_output.columns)[9:9+2]
fig, ax = plt.subplots()
x = np.arange(-40, 98, step=1)
print(x.shape)

for i in range(2):
    label = first_3[i]
    print(test_output[label])
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


# # Graph 5 Currents

# In[57]:


first_3 = list(test_output.columns)[11:11+3]
fig, ax = plt.subplots()
x = np.arange(-40, 98, step=1)
print(x.shape)

for i in range(3):
    label = first_3[i]
    print(test_output[label])
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


# # Graph 6 Voltages

# In[58]:


voltage = pd.DataFrame(pd.read_csv("../data_files/testV_output.csv"))
#print(voltage)
column = voltage.columns
print(voltage)
#print(voltage) # this shows the voltage 
first_3 =  voltage.columns
fig, ax = plt.subplots()
size = len(voltage[first_3[0]])
x = np.arange(0, size, step=1)


for i in range(1):
    label = first_3[i]
    print(label)
    #ax.scatter(x,test_output[label],color = color_list[i])
    ax.plot(x,voltage[label],color = color_list[i])
    
colorlist = zip(first_3,color_list)
handles = [mpatches.Patch(color=colour, label=first_3) for label, colour in colorlist]
t = ",".join(first_3)
ax.set_title(t)
plt.legend(handles, first_3, ncol=1, bbox_to_anchor=(1, 1))
plt.plot()
plt.tight_layout()
plt.savefig(f"../graphs/{t}.png")

