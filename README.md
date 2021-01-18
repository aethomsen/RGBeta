# RGBeta: Evaluator of Renormalization Group Beta-Functions
RGBeta is a Mathematica package that allows the user to extract RG beta-functions up to loop order 3-2-2 for gauge, Yukawa, and quartic couplings, respectively, for a large class of 4D renormalizable models.

Currently the project is in beta version. The accompanying paper and manual is forthcomming.  

## Installation
The simplest way to download and install RGBeta is to run the command
> Import["https://raw.githubusercontent.com/aethomsen/RGBeta/master/Install.m"]

to install RGBeta directly to the Applications folder in Mathematica's base directory, *$UserBaseDirectory*. After the install RGBeta can be loaded into any Mathematica notebook with
> << RGBeta`

As an alternative to the more permanent installation, simply download the github repository. RGBeta can then be run from a Mathematica notebook in the base directory with the following lines:
> SetDirectory@NotebookDirectory[];
> << RGBeta`


## Author
 - Anders Eller Thomsen (@aethomsen)

## License
RGBeta is free software under the MIT license.
