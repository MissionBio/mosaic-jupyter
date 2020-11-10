import os
import streamlit as st
import defaults as DFT

STATUS = None
ERROR = None
SUBHEADER = None


def rerun():
    raise st.script_runner.RerunException(st.script_request_queue.RerunData(None))


def subheader(msg):
    SUBHEADER.subheader(msg)


def info(msg, component=st, color=DFT.BLUE):
    component.markdown(f"<span style='color:{color}'>{msg}</span>", unsafe_allow_html=True)


def error(msg):
    global ERROR
    ERROR.markdown(f"<p style='font-size:18px'><span style='color:{DFT.RED}'><b>{msg}</b></span></p>", unsafe_allow_html=True)
    st.stop()


def status(msg):
    global STATUS
    STATUS.markdown(msg, unsafe_allow_html=True)


def download(download_path):
    name = download_path.split('/')[-1]
    os.popen(f'rm {DFT.DOWNLOADS_PATH}/*')
    os.popen(f'cp {download_path} {DFT.DOWNLOADS_PATH / name}')
    ERROR.markdown(f'<b>[Click here to download {name}](downloads/{name})</b>', unsafe_allow_html=True)


def init():
    global STATUS, ERROR, SUBHEADER

    st.title('Mosaic')
    ERROR = st.empty()
    SUBHEADER = st.empty()
    STATUS = st.empty()

    hide_streamlit_style = ('<style>\n'
                            '#MainMenu {visibility: hidden;}\n'
                            'footer {visibility: hidden;}\n'
                            '</style>')

    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

