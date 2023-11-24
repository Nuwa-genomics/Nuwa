import streamlit as st
import os
import subprocess

st.set_page_config(layout="wide", page_title='Nuwa', page_icon='🧬')

st.title("Terminal")

if 'wd' not in st.session_state:
    os.chdir(os.getenv('WORKDIR'))

with open('/app/css/common.css') as f:
    common_style = f"""
                <style>
                {f.read()}
                </style>
                """
    st.markdown(common_style, unsafe_allow_html=True)
    
with open('/app/css/terminal.css') as f:
    terminal_style = f"""
                <style>
                {f.read()}
                </style>
    """
    st.markdown(terminal_style, unsafe_allow_html=True)
    
bash, python = st.tabs(['Bash', 'Python'])

with bash:

    st.session_state['wd'] = "$"
    st.session_state['user'] = os.environ.get('HOME').split(sep='/')[-1]
    
    #subprocess.Popen(f"activate {os.getenv('CONDA_ENV')}", stdout=subprocess.PIPE, shell=True)

    col1, col2 = st.columns(2)
    container = st.container()
        
        
    cmd = col1.text_input(label="Bash terminal on host", key='ti_cmd', placeholder="Command")
                    
    if cmd:
        try:
            if cmd.find("cd") != -1:
                if len(cmd.split()) == 2:

                    st.session_state['stdout'] = f"""<div style='font: .9rem Inconsolata, monospace;'>{st.session_state.stdout}\n<div style='display: flex; gap:5px;'><div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$</div>{st.session_state.ti_cmd}</div></div>"""
                    os.chdir(cmd.split()[-1])
                    st.session_state['wd'] = os.getcwd()
                    container.markdown(f"""{st.session_state.stdout}\n<div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$<span class='cursor'>█</span></div>""", unsafe_allow_html=True)
            elif cmd.lstrip().lower() == "clear":
                st.session_state['stdout'] = f"""<div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$</div>"""
            else:
                bashCommand = st.session_state.ti_cmd
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE) #commands are passed through conda environment
                if 'stdout' not in st.session_state:
                    st.session_state['stdout'] = ""
                st.session_state['stdout'] = f"""<div style='font: .9rem Inconsolata, monospace;'>{st.session_state.stdout}\n{f"<div style='display: flex; gap:5px;'><div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$</div>{st.session_state.ti_cmd}</div>"}</div>"""
                container.markdown(st.session_state.stdout, unsafe_allow_html=True)
                full_output = ""
                output = ""
                while True:
                    with st.spinner("Executing command"):     
                        output = process.stdout.readline()

                        
                        # if stderr:
                        #     st.session_state['stdout'] = f"""<div style='font: .9rem Inconsolata, monospace;'>{st.session_state.stdout}\n{f"<div style='display: flex; gap:5px;'><div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$</div>{st.session_state.ti_cmd}</div>"}\n<div style='white-space: pre-line; color: red;'>{err.decode() for err in stderr}</div></div>"""
                        #     container.markdown(st.session_state.stdout, unsafe_allow_html=True)
                        #     break
                        
                        if not output:
                            break

                        output = output.decode()
                            
                        if 'stdout' not in st.session_state:
                            st.session_state['stdout'] = ""
                        container.markdown(f"""<div style='white-space: pre-line; font: .9rem Inconsolata, monospace;'>{output}</div></div>""", unsafe_allow_html=True)
                        full_output = f"{full_output}{output}"

                container.markdown(f"""<div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$<span class='cursor'>█</span></div>""", unsafe_allow_html=True)
                        
                st.session_state['stdout'] = f"""<div style='font: .9rem Inconsolata, monospace;'>{st.session_state.stdout}\n<div style='white-space: pre-line; font: .9rem Inconsolata, monospace;'>{full_output}</div></div></div>"""
                            
                        
        except Exception as e:
            if 'stdout' not in st.session_state:
                st.session_state['stdout'] = ""
            st.session_state['stdout'] = f"""<div style='font: .9rem Inconsolata, monospace;'>{st.session_state.stdout}\n{f"<div style='display: flex; gap:5px;'><div style='font-weight: bold; color: #1976d2;'>{st.session_state.user}@Nuwa:~{os.getcwd()}$</div>{st.session_state.ti_cmd}</div>"}\n<div style='white-space: pre-line; color: red;'>{e.__repr__()}</div></div>"""
            container.markdown(st.session_state.stdout, unsafe_allow_html=True)