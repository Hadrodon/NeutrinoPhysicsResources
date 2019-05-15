# NeutrinoPhysicsResources

Welcome\!

You're probably viewing this either through stumbling onto it accidentaly, viewing it intentionally, or your task is to continue my work.
If so, good luck, cause you are going to need it.

# Getting Started

To get started, you will need a version of ROOT installed to at least run the macros. The macros were written with ROOT v5.34.36, but more recent versions of ROOT can be used. To run a macro, you can simply type into the terminal `root (macroname)` or in the ROOT terminal, you can type `.x (macroname)` or `.L (macroname) $$ (functionname)`. This will run the macro, but some macros do need inputs to run.

A collection of files that can be used for the macros can be extracted from the .zip file for use in the macros. If you want to use your own macros, then you will need to first compile them using NUISANCE. You can build the CLAS6AcceptanceSummaryTree.so object using the shell script build.so.sh. You are going to need this shared object as well as NUISANCE to get the right object to compile the macros. This means that you will run NUISANCE with the file you want to process and the .so, along with an .xml with some acceptance measures to get this file. A sample of these xmls can be found in the CLAS6xmls folder for multiple files.

## Prerequisites

* You must be using a recent version of [CERN's ROOT software](https://root.cern.ch/). The macros were written with ROOT v5.34.36 but later versions can be used
* You must be using a recent version of [NUISANCE](https://nuisance.hepforge.org/)
* Given trees in zipped file are generated using the [GENIE-MC](http://www.genie-mc.org/) with version 2-12-2. If these are not used, you have to generate your own

## Debugging

* To debug any of the code, you will have to go into the macros yourself and change them manually to fix any bug found. If this is not possible or some things do not make sense, please consider contacting me at pontinicolas@gmail.com for further questions and debugging (Note: I will only be able to debug most of what I put up on GitHub and most likely won't be able to solve every problem you may encounter)

## Contributions

A big shoutout to Luke Pickering for his help in making the code. He gave the files necessary for the CLAS6AcceptanceSummaryTree as well as help for the physics behind event reconstruction. Without him there would be no CLAS6Acceptance. You can find his GitHub [here](https://github.com/luketpickering)