class component_base():
    def __init__(self, stylesheet):
        styles_path = "./components/styles/"
        with open(styles_path + stylesheet) as css:
            self.css = f"<style>{css.read()}</style>"

    def __repr__(self) -> str:
        return self.css + self.html