import streamlit as st
import pandas as pd

st.set_page_config(page_title='Nuwa', page_icon='🧬')

common_style = """
            <style>
            footer {visibility: hidden;}
            .st-emotion-cache-1cypcdb {background: linear-gradient(180deg, rgb(5, 39, 103) 0%, #3a0647 70%); box-shadow: 1px 0 10px -2px #000;}
            .st-emotion-cache-86cver {rgba(250, 250, 250, 0.6)}
            </style>
            """
st.markdown(common_style, unsafe_allow_html=True)


with open('css/workspace.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


st.title("Workspaces")


#conn = st.experimental_connection("postgresql", type="sql")

#df = conn.query('SELECT * FROM workspaces;', ttl="10m")
df = pd.read_csv('./sample_data.csv')

for row in df.iterrows():
    st.button(label=row[1][1], key=f"btn_{row[1][0]}")

#st.markdown(workspace(df), unsafe_allow_html=True)