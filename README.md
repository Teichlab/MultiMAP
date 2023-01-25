# MultiMAP
**MultiMAP** is a method for integrating single cell multi-omics. MultiMAP can also be used for batch correction. More detail is available in our [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02565-y).

<p align="center"><img src="docs/MultiMAP_schematic.png" width="900"></p>


## Installation

```bash
pip3 install git+https://github.com/Teichlab/MultiMAP.git
```

## Usage and Documentation

MultiMAP offers two functions accepting AnnData objects on input:
  - `MultiMAP.Integration()` expects a list of one AnnData per dataset, with the desired dimensionality reduction precomputed and stored in `.obsm`. This allows for refining the initial dimensionality reduction, e.g. if wishing to use `TFIDF_LSI` for ATAC data and PCA for RNA data.
  - `MultiMAP.Batch()` expects a single AnnData object with the dataset information stored in an `.obs` column. This allows for convenient integration with minimal preparation if all datasets can be treated with the same dimensionality reduction.

There's also an AnnData-independent `MultiMAP.matrix.MultiMAP()` function which operates directly on dimensionality reduction matrices. This requires precomputing all pairwise dimensionality reductions prior to calling MultiMAP.

A tutorial covering both RNA-ATAC integration and RNA-Seq batch correction use can be found [here](https://nbviewer.jupyter.org/github/Teichlab/MultiMAP/blob/master/examples/tutorial.ipynb).

Documentation of the function parameters can be found on [ReadTheDocs](https://multimap.readthedocs.io/en/latest/).

## Citation

If your work uses MultiMAP, please cite the [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02565-y):

	@article{jain2021multimap,
	  title={MultiMAP: dimensionality reduction and integration of multimodal data},
	  author={Jain, Mika Sarkin and Polanski, Krzysztof and Conde, Cecilia Dominguez and Chen, Xi and Park, Jongeun and Mamanova, Lira and Knights, Andrew and Botting, Rachel A and Stephenson, Emily and Haniffa, Muzlifah and others},
	  journal={Genome biology},
	  volume={22},
	  number={1},
	  pages={1--26},
	  year={2021},
	  publisher={BioMed Central}
	}

## Contact

Mika Sarkin Jain - mikasarkinjain@gmail.com \
Mirjana Efremova -  m.efremova@qmul.ac.uk \
Sarah Teichmann - st9@sanger.ac.uk
