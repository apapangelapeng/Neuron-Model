jupyter nbconvert --to python graph_test.ipynb #converts to python script
echo "converting display to print"
sed -i -e 's/display/print/g' graph_test.py
rm graph_test.py-e