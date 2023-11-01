# Description
C++ code Matthew used to generate the entrapment fasta files, which is an updated version of this repository from Lukas: https://github.com/percolator/mimic.

>With regard to the entrapment searches, here is a description about it from our picked group FDR manuscript:  
>For each entrapment database, we did an in-silico digest of the target database with trypsin-p as protease, without miscleavages and retained all peptides longer than 6 amino acids. We selected a fraction S of peptides to remain unchanged, thereby creating shared peptides between the target and entrapment databases. All other peptide sequences were shuffled, while keeping the C-terminal amino acid the same. Finally, these shuffled peptides replaced their original versions in the target protein sequence resulting in an entrapment protein, which then had the same number of shared peptides as the target version. We used two different ratios S of shared peptides, with S = 0.5 representing the shared ratio of the SwissProt+TrEMBL database and S = 0.04 representing the shared ratio of the SwissProt database.

> Some tips:
> - You probably don't need the part with the shared peptides at the moment, as this is to calibrate protein FDRs, so you could use S = 0.
>- I would suggest creating an entrapment database 9x the size of the target database. If your database is very large, 4x is probably also fine but might give less clear results in the comparison mentioned below.
>- For the search, you concatenate the target database with the entrapment database and use it as your new target database (in the mimic program, you can do this easily by using the -P flag). If you need to create a decoy database yourself, you will have to reverse the entire concatenated database (target+entrapment).
>- In the end you want to compare the decoy FDR (#decoys / #targets) with the entrapment FDR (#entrapments / #targets), which should be approximately equal.

# Conan
- install conan version 1.40 
- add new profile
```bash
conan profile new gcc --detect
conan profile update settings.compiler.libcxx=libstdc++11 gcc
```

# Building
To build the project please execute

```shell
$ cd src
$ cmake -DCMAKE_INSTALL_PREFIX=/tmp/install . && make && make install
```

# Run application
```shell
$ cd bin
$ ./mimic <path_to_fasta_file>
```
#### options

|Long flag         		|Short flag | description               												  			 |
|-----------------------|-----------|----------------------------------------------------------------------------------------|
|out-file          		|	 o      | File to write the mimicked protein sequences to (Default: prints to stdout) 			 |
|prefix			   		|    p      | Prefix to mimic proteins (Default: "mimic|Random_") 									 |
|mult-factor       		|    m      | sNumber of times the database should be multiplied (Default: 1) 						 |
|shared-pept-ratio      |    s      | Ratio of shared peptides that will stay preserved in the mimic database (Default: 0.0) |
|seed                   |    S      | Set seed of the random number generator. Default = 1									 |
|prepend                |    P      | Prepend the original fasta file to the output											 |
|help              		|    h      | produce help message      															 |


