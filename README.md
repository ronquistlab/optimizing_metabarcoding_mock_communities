# Evaluating metabarcoding protocols using replicated mock communities

Elzbieta Iwaszkiewicz-Eggebrecht, Emma Granqvist, Mateusz Buczek, Monika Prus, Tomas Roslin, Ayco J.M. Tack, Anders F. Andersson, Andreia Miraldo, Piotr Łukasik, Fredrik Ronquist

This is the source code repository accompanying the paper.

  * Preprint available at [bioRxiv](https://www.biorxiv.org/...)
  * Corresponding author: Fredrik Ronquist, fredrik.ronquist (at) nrm.se

## Contents

| Directory | Description                             |
|:----------|:----------------------------------------|
| `Birch`   | Models in Birch            |
| `ModelDescription`    | Detailed description of model and model extentions as PDF |
| `R`       | Help functions in R for data processing        |
| `WebPPL-exploratory`  | Early exploratory work with the models in WebPPL           |

## How to download this repository?

```
git clone http://github.com/ronquistlab/optimal_metabarcoding
```

# Models in Birch

## Birch Installation

[Birch](https://www.birch.sh) is run from the command line. Complete installation instructions are available on this [web page](https://www.birch.sh/getting-started/).

We have used the following Birch version: 1.634

Files in the [`config`](config) directory allow the configuration of various options such as the model, the number of samples to draw, and the number of particles to use. Files in the [`input`](input) directory provide the data sets in an appropriate format for use. Source files for all models can be found in the [`src`](src) directory.

## Running Birch code

Run an individual model with

```
birch sample --model OptmetabarModel_H --input input/H.json --output output/test-H.json --config config/configExp3-5.json
```

Output appears in the `output` directory, as specified in the command.

## Data input

In the [`input`](input) directory we provide the data sets in an appropriate format for use with Birch. "L" stands for lysate data and "H" for homogenate data. 
The datasets L_pred and H_pred contain the same data sets but in a slightly different format, as to suit the prediction code. 

## Different model versions

The basic model and its extensions is described in the [`ModelDescription`](ModelDescription) directory. 
The code either uses the lysate data set ("L") or homogenate data set ("H") and has been set up and used in four "versions". 

| Model name in Birch | Description                             |
|:----------|:----------------------------------------|
| `OptmetabarModel`   | The basic model, with one joint k and one joint theta for all spike-in species.            |
| `OptmetabarModel_4k`    | The model with one k per spike-in, and one joint theta for all spike-in species.  |
| `OptmetabarModel_4theta`       | The model with one joint k and one theta per spike-in species.       |
| `OptmetabarModel_4k4theta`  | The model with one k and one theta per spike-in species.          |
| `OptmetabarModel_4k4theta`  | The model with one k and one theta per spike-in species.          |

There are three folders in the  [`src`](src) directory, where the models in the N directory use number of spikein specimens, in the Weights directory the models use spikein weights, and in the Pred directory the code for predicting number of specimens or biomass is located. 


