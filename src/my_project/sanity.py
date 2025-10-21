import os
print("MPLBACKEND env:", os.environ.get("MPLBACKEND"))

import matplotlib
print("matplotlib backend:", matplotlib.get_backend())

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

x = np.linspace(0, 2*np.pi, 200)
y = np.sin(x)

fig, ax = plt.subplots()
ax.plot(x, y)
st.pyplot(fig)

st.success("Sanity plot rendered.")
