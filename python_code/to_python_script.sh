jupyter nbconvert --to python Static_AP.ipynb #converts to python script
echo "converting display to print"
sed -i -e 's/display/print/g' Static_AP.py
rm graph_test.py-e