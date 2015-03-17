### Overlapy
Project for the subjects of Structural Bioinformatics and Python (Master in _Bioinformatics for Health Sciences_ [Universitat Pompeu Fabra](http://www.upf.edu/bioinformatics/) - [Universitat de Barcelona](http://www.ub.edu/web/ub/ca/)); developed by Álvaro Abella Bascarán, Josep Arús Pous and Samuel Miravet Verde.

### What is Overlapy?

Overlapy is a python library and command line tool which provides a easy way to do protein superimposition. Overlapy is also provided in an even more user-friendly way, as a [web application](http://overlapy.undeadpixels.net).

If you are interested in know more about the program, its implementation and several examples of use, please read our [report](https://github.com/SMV818VMS/overlapy/blob/master/report/overlapy_report.pdf) (you have to download it as raw data in format pdb).

### Installation

Overlapy uses python3. Please check that it is installed:

```bash
python3 --version
```

If it's not, install it.

```bash
sudo apt-get install python3.4
```

In order to get Overlapy to work, you must first satisfy some library requirements. The preferred way is to use pip for this:

```bash
sudo apt-get install python3-pip
pip3 install -r requirements.txt
```

Once this is done, you can install the library + command line tool:

```bash
sudo python3 setup.py install
```

Now overlapy is available for you! You can use it as a python3 module, or as a command line tool:

#### Overlapy as a python module. Usage example
```python
import overlapy
# instantiate the Superimposer class
superimposer = overlapy.Superimposer()
# feed it the desired pdb files
superimposer.parse("file1.pdb", "file2.pdb")
# let's superimpose chains H and L with chains H and L
superimposer.select_chains(["H", "L"], ["H", "L"])
superimposer.superimpose()
# now let's get the results
rmsd = superimposer.rmsd
rotation_matrix = superimposer.rotation_matrix
alignment = superimposer.get_multiple_sequence_alignment()
# we can save the superimposed structure:
superimposer.save_superimposed_pdb("outputfile.pdb")
```

#### Overlapy as a command line tool.
Once you have run setup.py, you can execute it as a command line tool:

```bash
overlapy --help
```

To perform a basic superimposition just run:
```bash
overlapy -i1 protein1.pdb -i2 protein2.pdb -o protein1protein2.pdb
```

This will generate a PDB file containing the two initial structures superimposed (stored the input 1 as chain A and chain B for the second) and will show in your terminal the best Needleman & Wunsch alignment, the rotation matrix and the RMSD.

If you want to store in files the matrix (in tsv), the alignment (in clustal) and/or the RMSD (text file), just activate the options -r, -m and -a.

Finally, if you want to select chains to superimpose from the input PDBs, use the -c1 and -c2 arguments:

```bash
overlapy -i1 protein1.pdb -c1 HL -i2 protein2.pdb -c2 A -o protein1protein2.pdb -a -m
```
