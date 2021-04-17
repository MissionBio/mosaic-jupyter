import sys
from subprocess import check_call as sh

def replace_to_execute(nbname):

    finds = []
    replaces = []

    find = '"import missionbio.mosaic.io as mio\\n",'
    rep = '"import missionbio.mosaic.io as mio; '
    rep += 'import plotly.graph_objects as go; '
    rep += 'import plotly.offline as pyo; '
    rep += 'pyo.init_notebook_mode()\\n",'

    finds.append(find)
    replaces.append(rep)

    find = '"fig"'
    rep = '"go.Figure(fig)"'

    finds.append(find)
    replaces.append(rep)

    with open(nbname, 'r') as file:
        filedata = file.read()

    for find, replace in zip(finds, replaces):
        filedata = filedata.replace(find, replace)

    executable = nbname.replace('.ipynb', '.execute.ipynb')
    with open(executable, 'w') as file:
        file.write(filedata)

    return executable


def replace_back_in_notebook(nbname):

    finds = []
    replaces = []

    rep = '"import missionbio.mosaic.io as mio\\n",'
    find = '"import missionbio.mosaic.io as mio; '
    find += 'import plotly.graph_objects as go; '
    find += 'import plotly.offline as pyo; '
    find += 'pyo.init_notebook_mode()\\n",'

    finds.append(find)
    replaces.append(rep)

    rep = '"fig"'
    find = '"go.Figure(fig)"'

    finds.append(find)
    replaces.append(rep)

    with open(nbname, 'r') as file:
        filedata = file.read()

    for find, replace in zip(finds, replaces):
        filedata = filedata.replace(find, replace)

    with open(nbname, 'w') as file:
        file.write(filedata)


def convert_nb(nbname):

    executable = replace_to_execute(nbname)

    # Execute the notebook
    sh(["jupyter", "nbconvert", "--to", "notebook",
        "--execute", "--inplace", executable])

    replace_back_in_notebook(executable)

    # Convert to .html
    output = nbname.replace('.ipynb', '')
    sh(["jupyter", "nbconvert", "--to", "html_ch", executable, "--output", output])
    sh(["rm", executable])


if __name__ == "__main__":

    for nbname in sys.argv[1:]:
        convert_nb(nbname)
