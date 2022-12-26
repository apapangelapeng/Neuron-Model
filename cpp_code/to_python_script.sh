jupyter nbconvert --to python ../python_code/Static_AP.ipynb #converts to python script
echo "converting display to print"
sed -i -e 's/display/print/g' ../python_code/Static_AP.py
rm  ../python_code/Static_AP.py-e
python3 ../python_code/Static_AP.py