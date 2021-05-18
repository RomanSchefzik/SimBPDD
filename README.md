# SimBPDD

Simulating differential distributions for Beta-Poisson models, in particular for single-cell RNA-sequencing (scRNA-seq) data

> For more details on the methods and further information, please have a look at the accompanying paper:  
> R. Schefzik (2021). SimBPDD: Simulating differential distributions in Beta-Poisson models, in particular for single-cell RNA-sequencing data. *Annales Mathematicae et Informaticae*, 53:283-298. Available at https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_53_from283to298.pdf

## Installation

### Prerequisites

Please make sure that you have the latest version of R installed. To use the SimBPDD package, it is required to install the following packages first:  
`install.packages("statmod")`  
`install.packages("doParallel")`  
`install.packages("foreach")`  

### Installation of the SimBPDD package

To install the SimBPDD package, run the following:  
`install.packages("devtools")`  
`devtools::install_github("RomanSchefzik/SimBPDD")` 

The SimBPDD package can then be loaded using  
`library(SimBPDD)`

## Usage

For the usage of the package and a detailed description of the functions, see the accompanying manual.
