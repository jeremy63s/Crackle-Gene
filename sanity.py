import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

st.write("matplotlib backend:", matplotlib.get_backend())

x = np.linspace(0, 2*np.pi, 200)
fig, ax = plt.subplots()
ax.plot(x, np.sin(x))
st.pyplot(fig)
st.success("Sanity OK")
