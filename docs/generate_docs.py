import ast
import os
import glob
import re
from docstring_parser import parse

########################################################################################################
# This is a script to generate documentation from docstrings. Must be run from 'Nuwa' directory.       #
# NOTE: Automatically run as part of the github actions workflow and does not need to be run manually. #
########################################################################################################

src_root = "app/pages" #root folder for source code

python_files = glob.glob("*.py", root_dir=src_root)

for file in python_files:
        with open(os.path.join(src_root, file)) as fd:
            print()
            print(file)
            print("Functions:")
            file_contents = fd.read()
            module = ast.parse(file_contents)
                
            function_definitions = [node for node in module.body if isinstance(node, ast.FunctionDef)]
            for f in function_definitions:
                print(f.name)
                if ast.get_docstring(f):
                    print(ast.get_docstring(f))
                else:
                     print("WARNING: No Doctype")
                

            class_definitions = [node for node in module.body if isinstance(node, ast.ClassDef)]
            for class_ in class_definitions:
                print("\nClass name:", class_.name)
                methods = [n for n in class_.body if isinstance(n, ast.FunctionDef)]
                for method in methods:
                    print(method.name)
                    if ast.get_docstring(method):
                        docstring_raw = ast.get_docstring(method)
                        docstring = parse(docstring_raw)
                        markdown = f"# {method.name}\n{docstring.short_description}"

                        #parameters
                        if len(docstring.params) > 0:
                            markdown = markdown + "\n## Parameters"
                            for i, param in enumerate(docstring.params):
                                markdown += f"\n```{docstring.params[i].arg_name}: {docstring.params[i].type_name}``` {docstring.params[i].description} \n<hr style='height:1px'>"

                        #screenshots
                        if len(docstring.examples) > 0:
                            markdown += "\n## Web view"
                            for i, example in enumerate(docstring.meta):
                                #look in notes heading
                                if docstring.meta[i].args[0] == 'notes':    
                                    #capture image string according to numpy docstring format
                                    image_md = re.split(r'image::\s*', docstring.meta[i].description)[-1]
                                    markdown += f"\n<img alt='{method.name}_screenshot' src='{image_md}' width='300' height='500'>"
                    
                        
                        #python equivalent
                        markdown = markdown + "\n## Python equivalent"
                        for i, example in enumerate(docstring.examples):
                            markdown = markdown + f"\n```python\n{docstring.examples[i].description}\n```"
                        
                        
                        f = open(f"./docs/reference/{class_.name}/{method.name}.md", "w")
                        f.write(markdown)
                        f.close()
                    else:
                        print("WARNING: No Docstring")

                    

