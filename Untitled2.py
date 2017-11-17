
# coding: utf-8

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth


# In[7]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

fp = 'C:\Users\sleblan2\Research\ORACLES\data_other\\'


# In[4]:

coma = load_ict(fp+'COMA_P3_20160920_R0.ict')


# In[5]:

cas = load_ict(fp+'UND-CAS_P3_20160920_R1.ict')


# In[20]:

p3 = load_ict(fp+'Hskping_P3_20160920_R1.ict')


# In[21]:

plt.figure()
plt.plot(cas['Start_UTC'],cas['CAS_NtAero'])
plt.plot(coma['Start_UTC'],coma['CO_ppbv'],'g')
plt.plot(p3['Start_UTC'],p3['Static_Pressure']/10,'r')


# In[17]:

0.006*60*60


# In[19]:

(9.64211-9.63748)*3600


# In[22]:

(8.23844-8.23477)*3600


# In[23]:

(8.74984-8.74585)*3600


# In[24]:

0.0068*3600


# In[ ]:



