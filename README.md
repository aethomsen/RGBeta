# RGBeta: Evaluator of Renormalization Group Beta-Functions
RGBeta is a Mathematica package that allows the user to extract RG beta-functions up to loop order 3-2-2 for gauge, Yukawa, and quartic couplings, respectively, for a large class of 4D renormalizable models. 

Currently the project is in beta version. The accompanying paper and manual is forthcomming.  

## Installation 
RGBeta can be run from a Mathematica notebook in the base directory with the following lines
> SetDirectory@NotebookDirectory[];
> << RGBeta`

Alternatively, one can load it from anywhere by providing the path to the directory of RGBeta with   
> AppendTo[$Path, " Directory "];
> << RGBeta`

using the appropriate directory.


## Author
 - Anders Eller Thomsen (@aethomsen) 

## License
RGBeta is free software under the MIT license. 
