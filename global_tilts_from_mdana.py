#!/usr/bin/env python
# coding: utf-8

# In[4]:


import MDAnalysis as mda
from MDAnalysis.analysis import helix_analysis as hel
import pandas as pd
import nglview as nv
from IPython.display import display
import warnings
import os

get_ipython().run_line_magic('matplotlib', 'inline')


# In[5]:


# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")


# In[11]:



pdb_files = ["8et6.pdb", 
             "8jts.pdb", 
             "8jtt.pdb", 
             "8jtv.pdb", 
             "8sc4.pdb", 
             "8sc1.pdb"]


# In[12]:


# Define TM helix ranges (residue numbers)
tm_ranges = [
    (15, 42),    # TM1
    (146, 170),  # TM2
    (177, 195),  # TM3
    (200, 228),  # TM4
    (233, 257),  # TM5
    (261, 283),  # TM6
    (341, 368),  # TM7
    (377, 397),  # TM8
    (403, 423),  # TM9
    (431, 455),  # TM10
    (462, 490),  # TM11
    (493, 513)   # TM12
]


# In[13]:


def analyze_helix_tilts(inp_pdb, tm_ranges):
    """Analyze tilt angles for all transmembrane helices"""
    u = mda.Universe(inp_pdb)
    results = {'system': os.path.basename(inp_pdb)[8:-4]}
    
    for idx, (start, end) in enumerate(tm_ranges, 1):
        try:
            h = hel.HELANAL(
                u,
                select=f'name CA and resnum {start}-{end}',
                ref_axis=[0, 0, 1]
            ).run()
            
            # Store rounded results
            tilt_mean = round(h.results.summary['global_tilts']['mean'], 2)
            results[f'tm{idx}_tilt'] = tilt_mean
            
        except Exception as e:
            print(f"Error processing TM{idx} ({start}-{end}): {str(e)}")
            results[f'tm{idx}_tilt'] = None
    
    return results


# In[15]:


# Create DataFrame
df = pd.DataFrame([analyze_helix_tilts(f, tm_ranges) for f in pdb_files])

# # Save to CSV
df.to_csv('oct1_tm_helix_tilts.csv', index=False)

# Display results
print("\nFinal DataFrame:")
display(df)


# In[ ]:




