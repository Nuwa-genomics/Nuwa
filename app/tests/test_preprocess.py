from streamlit.testing.v1 import AppTest
import time

class Test_Preprocess:
    def __init__(self):
        at = AppTest.from_file("Dashboard.py")
        at.run()
        assert not at.exception