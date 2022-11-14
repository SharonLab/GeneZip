# GeneZip

This program builds LZ78 models for DNA sequences from a training set that consists of fasta files and calculates the expected log-loss for a testing set of (other) fasta files.

## Installation

### Installation option 1: Using a prebuilt binary

Linux users can download a prebuilt release from [here](https://github.com/SharonLab/GeneZip/releases/latest/download/GeneZip.tar.gz).
The prebuilt binary was tested on Ubuntu 22.04.

Run the following commands from the directory in which you want GeneZip to be installed:

```bash
wget https://github.com/SharonLab/GeneZip/releases/latest/download/GeneZip.tar.gz
tar -xf GeneZip.tar.gz
cd GeneZip
./GeneZip --help
```
To make GeneZip available from anywhere in your system, add the GeneZip binary directory path to the PATH environment variable.

### Installation option 2: installing from source

 1. Clone the repository using `git clone git@github.com:SharonLab/GeneZip.git`.
 2. Change directory into the code directory using `cd GeneZip/code`.
 3. Either use pre-installed rust's cargo to build the project with `cargo build --release`, or use the provided build.sh script to install rust+cargo and build the project with `./build.sh`.
 4. The binary file will be located within a directory called target/release and can be executed using `./target/release/GeneZip --help`.
 5. To access GeneZip from anywhere in your system, copy the file ./target/release/GeneZip into a folder that is under the environment variable PATH, or add GeneZip's directory to the PATH variable.



## GeneZip examples

Two small examples are included in the "examples" directory. 
Change your working directory to the examples/tiny_example or to examples/small_example, and run the run.sh file. 

```bash
./run.sh
```

The output files will be in the 'output' folder.

## Running the program
We recommend that the program will be used with gene files for training and testing rather than genome files because gene files achieve better classification results. For genomes deposited in NCBI you can usually download also gene files. You can also use gene prediction programs such as prodigal to create gene files from assembled genomes. 

Program command line is
```
GeneZip -i <training-name2file> -t <predict-name2file> -o <output> [-d <max-depth>] [-j <max-jobs>]
```
Where
- **-i   \<training-name2file\>** The training set file. It consists of two tab separated columns:<br>
  \<cluster-name\>  \<path-to-fasta-file\><br>
  The fasta file may contain one or more sequences.
  Multiple lines may have the same cluster name, in which case all the files will be used to train the model for the same cluster.<br>
- **-t   \<predict-name2file\>** The testing set. It consists of two tab separated columns:<br>
  \<ID\> \<path-to-fasta-file\><br>
  The ID will be used for identification of the tested fasta in the output. Identical IDs will not be clustered. The fasta file may contain a genome, some genes, etc.
- **-o  \<output\>** Output file name
- **-d  \<depth\>** The maximum depth for the LZ78 tree. Less than 9 will result with poor results, more than 13 may require too much disk space. This flag is optional, the default depth is set to 13 and should work well if one genome per cluster is used.
- **-j  \<max-jobs\>** The upper limit of threads to use. GeneZip makes a good use of multithreading so consider allowing it to use as many threads as possible, especially if you are working on large training or testing sets (more than 20-30 fasta files each). Uses 1 thread by default, set to 0 to allow GeneZip to use all available threads.


Two example files and commands are provided in the [examples/](/examples/) directory

In both examples, the script 'run.sh' will run the example from within the same folder. The files 'training.txt' and 'testing.txt' represent \<training-name2file\> and \<predict-name2file\> (respectevly).

The data used in the examples can be found in the [data/](/data/) folder. All included files are fasta files.

## Results
The program writes a table to \<output\> . For each file in the test set (row), the log-loss of all models in the training set (columns) are written. The program also writes statistics about the different models to the standard error, this can be piped to a file using `2> filename`.


## C version

The folder c_version contains the unsupported version of the project developed in C.
It contains a make file and two small running examples.

## How to cite and ask for support
The GeneZip project was presented in the [ISBRA 2022](https://mangul-lab-usc.github.io/ISBRA/) conference and will be published as a paper soon.

For support, please contact Yochai Meir at yochai AT titat DOT info or open an issue.

