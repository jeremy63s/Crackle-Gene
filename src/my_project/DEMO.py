# UI.py
import io
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Streamlit Demo Playground", layout="wide")

# ---------------- Sidebar ----------------
with st.sidebar:
    st.header("Controls")
    name = st.text_input("Your name", value="Guest")
    slider_n = st.slider("Rows for demo DataFrame", 10, 1000, 200, step=10)
    noise = st.slider("Noise level", 0.0, 1.0, 0.2, step=0.05)
    st.caption("Change values and watch the app rerun.")

# ---------------- Session State ----------------
if "clicks" not in st.session_state:
    st.session_state.clicks = 0

def bump():
    st.session_state.clicks += 1

st.title("YOU ARE RUNNING DEMO")
st.write(f"Hello, **{name}**.")

st.button("Increment counter", on_click=bump)
st.write(f"Counter value: `{st.session_state.clicks}`")

st.divider()

# ---------------- Tabs ----------------
tab1, tab2, tab3, tab4 = st.tabs(["Widgets + Form", "Data & Charts", "Files", "State & Caching"])

# ---- Tab 1: Widgets + Form ----
with tab1:
    st.subheader("Widgets & a Form")
    col1, col2, col3 = st.columns(3)
    with col1:
        choice = st.selectbox("Pick an option", ["Alpha", "Beta", "Gamma"])
    with col2:
        agree = st.checkbox("I agree")
    with col3:
        level = st.radio("Level", ["Low", "Medium", "High"], index=1)

    with st.form("feedback_form", clear_on_submit=True):
        st.write("Quick feedback form")
        mood = st.text_input("How are you feeling?")
        submitted = st.form_submit_button("Submit")
        if submitted:
            st.success(f"Thanks, recorded: {mood!r}")

    with st.expander("See code of this section"):
        st.code("""
choice = st.selectbox(...)
agree = st.checkbox(...)
level = st.radio(...)
with st.form("feedback_form"):
    mood = st.text_input("How are you feeling?")
    st.form_submit_button("Submit")
        """, language="python")

# ---- Tab 2: Data & Charts ----
with tab2:
    st.subheader("Data & Charts")
    x = np.linspace(0, 4*np.pi, slider_n)
    y = np.sin(x) + noise * np.random.randn(slider_n)
    df = pd.DataFrame({"x": x, "y": y})
    st.dataframe(df.head(10), width=True)

    st.write("Line chart (built-in):")
    st.line_chart(df.set_index("x"))

    st.write("Matplotlib scatter (custom):")
    fig, ax = plt.subplots()
    ax.scatter(df["x"], df["y"])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Scatter of noisy sin(x)")
    st.pyplot(fig)

    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button("Download CSV", data=csv_bytes, file_name="demo_data.csv", mime="text/csv")

# ---- Tab 3: Files ----
with tab3:
    st.subheader("File upload & preview")
    uploaded = st.file_uploader("Upload a CSV or TXT file", type=["csv", "txt"])
    if uploaded:
        if uploaded.type == "text/csv" or uploaded.name.lower().endswith(".csv"):
            try:
                df_up = pd.read_csv(uploaded)
                st.write("Detected CSV. Preview:")
                st.dataframe(df_up.head(), width=True)
            except Exception as e:
                st.error("Could not read as CSV.")
                st.exception(e)
        else:
            text = uploaded.getvalue().decode("utf-8", errors="replace")
            st.write("Detected text file. First 400 chars:")
            st.code(text[:400])

# ---- Tab 4: State & Caching ----
with tab4:
    st.subheader("Caching demo")

    @st.cache_data(show_spinner=False)
    def slow_square(n: int) -> int:
        time.sleep(0.5)  # pretend this is expensive
        return n * n

    num = st.number_input("Enter an integer", min_value=0, max_value=10000, value=12, step=1)
    st.write("Result (cached):", slow_square(int(num)))

    st.info("Try changing the value—cached results are instant when inputs are the same.")

st.caption("Tip: use st.session_state to persist values; use st.cache_data for pure functions and st.cache_resource for heavy stateful objects.")
