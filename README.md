Notas acerca de prepare_aromatic_center.py

Este programa tiene la capacidad de devolver uno o más archivos pdbqts junto con los nuevos centros de masa de cada anillo aromático. Antes de empezar deberá tener instalado open babel, pybel y numpy.

El modo de ejecución es el siguiente:
 $ ./prepare_aromatic_center.py -l loratadine.pdbqt -o loratadine_dummy.pdbqt
si no especificamos la salida por defecto escribe "filename" + _dummy.pdbqt. Por otra parte si el archivo de entrada es un sdf, toma el nombre de cada molécula de salida de la sección "title", y el nombre de salida es similar al descripto anteriormente.

El script usted lo puede obtener File:Prepare aromatic center.zip aqui junto con un sdf de prueba. Al ejecutar el archivo de prueba obtendrá:

 $ ./prepare_aromatic_center.py -l SMILES.sdf 
 WARNING: the input file  azulene_dummy.pdbqt  does not contain aromatic rings
 WARNING: the input file  ciclopropene_dummy.pdbqt  does not contain aromatic rings

Problemas reportados:

El uso de anillos aromáticos fusionados puede llevar a que anillos individuales no sean aromáticos, pero el sistema completo si lo sea. Un ejemplo de esto es azuleno.

Picture 8.png

Para este sistema el programa no retorna centros de masa (CM) para los anillos, puesto que el algoritmo GetSSSR del open babel identifica cada uno de ellos como No aromático
# prepare_aromatic_center
