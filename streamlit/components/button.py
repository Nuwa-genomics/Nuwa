from components.component_base import component_base

class button(component_base):   
    def __init__(self, text, icon=None, type='solid', width='100%'):
        class_name = 'button'
        super().__init__(stylesheet='button.css')
        self.html = f" \
        <button style='width: {width}' class={class_name}> \
            {text} \
        </button>"