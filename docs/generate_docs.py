import ast
import os
import glob
from docstring_parser import parse

root = "app/pages"

python_files = glob.glob("*.py", root_dir=root)
print(python_files)



for file in python_files:

        with open(os.path.join(root, file)) as fd:
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
                        markdown = f"## {method.name}"
                        docstring_raw = ast.get_docstring(method)
                        docstring = parse(docstring_raw)

                        #parameters
                        if len(docstring.params) > 0:
                            markdown = markdown + "## Parameters"
                            for i, param in enumerate(docstring.params):
                                markdown += f"```{docstring.params[i].arg_name}: {docstring.params[i].type_name}```"
                                markdown += docstring.params[i].description
                        
                        f = open(f"./docs/reference/{class_.name}/{method.name}.md", "w")
                        f.write(f"## {method.name}\n{ast.get_docstring(method)}")
                        f.close()
                    else:
                        print("WARNING: No Doctype")

            print()
            print()
            print()
                    

