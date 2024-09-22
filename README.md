# Scripts para Validación Cruzada y HMMsearch en Lote

Este repositorio contiene scripts diseñados para la automatización del procesamiento de secuencias biológicas con el software HMMER. 
El script `k_fold_cross_validation.py` permite la creación de modelos ocultos de Markov (HMM) y realización de validaciones cruzadas de k iteraciones a un conjunto de secuencias de interés.
El script `batch_hmmsearch.py` ejecuta en lote un escaneo de un HMM a un conjunto de archivos .FASTA, aislado las secuencias encontradas y cuantificándolas por archivo.

## Requisitos del Sistema

- Sistema operativo: Linux
- Dependencias de software externo:
  - [hmmer](http://hmmer.org/) (verificado en versión 3.3.2)
  - [clustalo](http://www.clustal.org/omega/) (verificado en versión 1.2.4)
  - [seqkit](https://bioinf.shenwei.me/seqkit/) (verificado en versión 2.1.0)
  - [cd-hit](cd-hit (https://github.com/weizhongli/cdhit/tree/master)) (verificado en versión 4.8.1)
- Paquetes de Python:
  - [Biopython](https://biopython.org/) (verificado en versión 1.84)

**Nota**: 
1. Instalar `hmmer`, `clustalo` y `seqkit` de manera separada y asegurarse de que estén accesibles en el PATH del sistema. También asegúrate de tener `Biopython` instalado en tu entorno de Python.
2. La clase `HmmDataProcessor`, contenida en este repositorio, es requerida para la ejecución de los scripts debido a que posee funciones para el procesamiento del resultado del hmmsearch que es utilizado por ambos scripts.

## Scripts

### 1. Validación Cruzada de k iteraciones

Este script ejecuta una validación cruzada de k iteraciones para datos de secuencias de biológicas con diferentes niveles de redundancia, que se evalúan como un hiperparámetro. También requiere un conjunto de datos negativos para pruebas.
Para cada redundancia establecida retorna una carpeta que contiene la data posterior a la reducción de redundancia, la partición de la secuencias para las k iteraciones, los alineamientos múltiples de secuencia, los HMM y los resultados del `hmmsearch` en cada iteración.
Adicionalmente crea un archivo `results.csv` que cuantifica los verdaderos positivos (VP), falsos negativos (FN), verdaderos negativos (VN) y falsos positivos (FP) de todos los modelos creados.

#### Uso

```bash
python3 k_fold_cross_validation.py <nombre_trabajo> -ts <archivo_entrenamiento.fasta> -ns <archivo_negativo.fasta> -k <k_folds> -r <redundancia>
```

#### Ejemplo

```bash
python3 k_fold_cross_validation.py GAF_domain -ts GAF_training_data.fasta -ns negative_data.fasta -k 5 -r 100,90,80,70,60
```

#### Parámetros
- job_name: Nombre del trabajo.
- -ts, --training_data_file: Archivo de datos de entrenamiento en formato FASTA.
- -ns, --negative_set_file: Archivo de conjunto negativo en formato FASTA.
- -k, --k_folds: Número de particiones para la validación cruzada (predeterminado: 5).
- -r, --redundancy: Niveles de redundancia a evaluar como hiperparámetro, separados por comas (predeterminado: 100).

### 2. Búsqueda HMM en Lote

Este script ejecuta `hmmsearch` con un HMM dado en un lote de archivos FASTA especificados en un archivo CSV.
Cada archivo analizado genera una carpeta que contiene el resultado del `hmmsearch` y un archivo FASTA con las secuencias encontradas bajo el umbral de E-value establecido.
Además, crea un archivo `quantification.csv` con la cantidad de proteínas encontradas en cada archivo FASTA.

#### Uso

```bash
python3 batch_hmmsearch.py <nombre_trabajo> <HMM.hmm> <archivos_fasta.csv> <valor_e>
```

#### Ejemplo

```bash
python3 batch_hmmsearch.py GAF HMM.hmm FASTA_files.csv 1e-5
```

#### Parámetros
- job_name: Nombre del trabajo.
- hmm_path: Ruta al archivo HMM.
- FASTA_files: Archivo CSV con las rutas de los archivos FASTA.
- E_value: Umbral de valor E para las búsquedas.
