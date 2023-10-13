from components.component_base import component_base
from components.workspace_tile import workspace_tile

class workspace(component_base):   
    def __init__(self, df):
        class_name = 'workspace'
        super().__init__(stylesheet='workspace.css')
        self.html = ""

        # tiles is an array containing: name, icon, date created
        for tile in df.iterrows():
            self.html += workspace_tile(name=tile[1][1], icon='dna').html

        self.html += workspace_tile(name='New workspace', icon='plus').html

        self.html = f"<div class={class_name}>{self.html}</div>"