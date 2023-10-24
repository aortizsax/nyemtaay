# Nyemtaay

## Introduction

Nyemtaay is a Python library for information theoretic gene flow analysis and population genetic calculations.
The name ``nyemtaay`` is from a [Kumeyaay](https://en.wikipedia.org/wiki/Kumeyaay_language) ([Ipai](https://en.wikipedia.org/wiki/Ipai_language)) word for "[cougar](https://livingdictionaries.app/iipay-aa/entry/TF63My3zWPXP1qElSHcH)", in homage to the indigenous peole of Southern California, on whose land we live and work.

## Requirements and Installation
Nyemtaay runs under python3 (>3.9)

We recommend that you install directly from the main GitHub repository using pip (which works with an Anaconda environment as well):

```
$ python3 -m pip install --user --upgrade git+https://github.com/aortizsax/nyemtaay.git
```

## Applications

### Analysis
Once cloned locally, the trial data can be check using

	$ cd nyemtaay/
	$ nyemtaay -m ./tests/testdata/*.csv -f ./tests/testdata/*.fasta
