import os
import sys
from streamlit import cli as stcli

if __name__ == '__main__':
    launchdir = os.path.dirname(sys.argv[0])
    sys.argv = ["streamlit", "run", f"{launchdir}/mosaic.py", "--global.developmentMode=false"]
    sys.exit(stcli.main())
