from components.component_base import component_base

class workspace_tile(component_base):   
    def __init__(self, name, icon, type='solid'):
        class_name = 'workspace_tile'
        super().__init__(stylesheet='workspace.css')
        self.html = f" \
        <link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css'> \
        <div class={class_name}> \
            {name} \
            <i class='fa-{type} fa-{icon}'></i> \
        </div>"
    