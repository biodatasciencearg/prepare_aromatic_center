Notas acerca de prepare_aromatic_center.py

Este programa tiene la capacidad de devolver uno o más archivos pdbqts junto con los nuevos centros de masa de cada anillo aromático. Antes de empezar deberá tener instalado open babel, pybel y numpy.

para instalar pybel (ubuntu/debian)

apt-get install python-openbabel

SE RECOMIENDA FUERTEMENTE EL EMPLEO DE ARCHIVOS SDF QUE CUENTAN CON MUCHAS MAS INFORMACION QUE LOS ARCHIVOS PDQBTs.

El modo de ejecución es el siguiente:
 $ ./prepare_aromatic_center.py -l loratadine.pdbqt -o loratadine_dummy.pdbqt
si no especificamos la salida por defecto escribe "filename" + _dummy.pdbqt. Por otra parte si el archivo de entrada es un sdf, toma el nombre de cada molécula de salida de la sección "title", y el nombre de salida es similar al descripto anteriormente.

El script usted lo puede obtener File:Prepare aromatic center.zip aqui junto con un sdf de prueba. Al ejecutar el archivo de prueba obtendrá:

 $ ./prepare_aromatic_center.py -l SMILES.sdf 
 WARNING: the input file  azulene_dummy.pdbqt  does not contain aromatic rings
 WARNING: the input file  ciclopropene_dummy.pdbqt  does not contain aromatic rings

