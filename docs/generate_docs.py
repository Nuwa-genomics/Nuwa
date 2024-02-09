import ast
import os
import glob
import re
from docstring_parser import parse
import shutil

########################################################################################################
# This is a script to generate documentation from docstrings. Must be run from 'Nuwa' directory.       #
# NOTE: Automatically run as part of the github actions workflow and does not need to be run manually. #
########################################################################################################

if os.path.exists('./docs/reference'):
    shutil.rmtree('./docs/reference') #remove old

# create root folder for API reference
os.mkdir('./docs/reference')

#write main mardown file
with open(f"./docs/reference/README.md", "w") as f:
    f.write("# Reference\n\n{% include list.liquid all=true %}")

classes_ordered = ["Dashboard", "Upload", "Preprocess", "Integrate", "Create_CiteSeq_model", "Create_Solo_model", "Cluster_analysis", "Differential_gene_expression", "Trajectory_inference", "Spatial_transcriptomics", "Terminal", "Plotly_3D", "Utils"]

#make the sub dirs
for cls in classes_ordered:
    os.mkdir(f"./docs/reference/{cls}")


src_root = "app"
page_root = "app/pages" #root folder for source code

python_files = glob.glob("*.py", root_dir=page_root, recursive=False)
dashboard_file = glob.glob("*.py", root_dir=src_root, recursive=False)
python_files = dashboard_file + python_files

for file in python_files:
        path = os.path.join(src_root, file) if file == "Dashboard" else os.path.join(page_root, file)
        with open(path) as fd:
            file_contents = fd.read()
            module = ast.parse(file_contents)
                
            #functions not contained in a class
            function_definitions = [node for node in module.body if isinstance(node, ast.FunctionDef)]
            for f in function_definitions:
                if ast.get_docstring(f):
                    #print(ast.get_docstring(f))
                    raise NotImplemented
                else:
                     print(f"WARNING: No Doctype: {f.name} in {file}")
                

            class_definitions = [node for node in module.body if isinstance(node, ast.ClassDef)]
            for class_ in class_definitions:

                if os.path.exists(f'./docs/reference/{class_.name}'):
                    with open(f"./docs/reference/{class_.name}/README.md", "w") as f:
                        markdown = f"---\nsort: {classes_ordered.index(class_.name)}\n---\n"
                        markdown += f"# {class_.name.replace('_', ' ')}\n\n"
                        markdown += f"{parse(ast.get_docstring(class_)).short_description}\n\n"
                        #screenshots
                        docstring = parse(ast.get_docstring(class_))
                        for i, example in enumerate(docstring.meta):
                            #look in notes heading
                            if docstring.meta[i].args[0] == 'notes':    
                                #capture image string according to numpy docstring format
                                images = re.split(r'..\s*image::\s*', docstring.meta[i].description)
                                for img in images:
                                    if img.__contains__("screenshot"):
                                        markdown += f"\n<img style='border-radius:15px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='{img.strip()}'>"
                        #methods
                        markdown += f"\n## Methods:\n\n{{% include list.liquid all=true %}}"
                        f.write(markdown)

                methods = [n for n in class_.body if isinstance(n, ast.FunctionDef)]
                for method in methods:
                    if ast.get_docstring(method):
                        docstring_raw = ast.get_docstring(method)
                        docstring = parse(docstring_raw)
                        markdown = f"# {method.name.replace('_', ' ')}\n{docstring.short_description}"

                        #parameters
                        if len(docstring.params) > 0:
                            markdown = markdown + "\n## Parameters"
                            for i, param in enumerate(docstring.params):
                                markdown += f"\n```{docstring.params[i].arg_name}: {docstring.params[i].type_name}``` {docstring.params[i].description} \n<hr style='height:1px'>"

                        #screenshots
                        markdown += "\n## Web view"
                        for i, example in enumerate(docstring.meta):
                            #look in notes heading
                            if docstring.meta[i].args[0] == 'notes':    
                                #capture image string according to numpy docstring format
                                images = re.split(r'..\s*image::\s*', docstring.meta[i].description)
                                for img in images:
                                    if img.__contains__("screenshot"):
                                        markdown += f"\n<img style='border-radius:15px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='{method.name}_screenshot' src='{img.strip()}'>"
                    
                        
                        #python equivalent
                        if len(docstring.examples) > 0:
                            markdown = markdown + "\n## Python equivalent"
                            for i, example in enumerate(docstring.examples):
                                markdown = markdown + f"\n```python\n{docstring.examples[i].description}\n```"
                        
                        
                        f = open(f"./docs/reference/{class_.name}/{method.name}.md", "w")
                        f.write(markdown)
                        f.close()
                    else:
                        print(f"WARNING: No Docstring: {class_.name}.{method.name}")

                    

