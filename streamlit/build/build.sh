pyinstaller --additional-hooks-dir=hooks \
			--hidden-import missionbio.mosaic \
			--hidden-import missionbio.h5 \
			--hidden-import plotly \
			--clean --onefile run_mosaic.py

mv ./dist/* ../
rm -r ./build ./run_mosaic.spec ./dist
