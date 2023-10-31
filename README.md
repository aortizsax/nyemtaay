# Nyemtaay

## Introduction

Nyemtaay is a python library for population genetics caluclations. 
The name ``nyemtaay`` is from a [Kumeyaay](https://en.wikipedia.org/wiki/Kumeyaay_language) ([Ipai](https://en.wikipedia.org/wiki/Ipai_language)) word for "[cougar](https://livingdictionaries.app/iipay-aa/entry/TF63My3zWPXP1qElSHcH)", in homage to the inigenous peole of Southern California, on whose land we live and wok. 

## Requirements and Installation
Nyemtaay runs under python3 (>3.9) 

We recommend that you install directly from the main GitHub repository using pip (which works with an Anaconda environment as well):

```
$ python3 -m pip install --user --upgrade git+https://github.com/aortizsax/nyemtaay.git
```

## Applications

Below I will outline the applications icluding traditional popgen stats, information theory analysis, and tissue of origin for a cancer dataset. 

Once cloned locally, check the help page

```
nyemtaay -h
```

### Information Theory   

The greek bears data will be used to show the functionality of Nyemtaay. 

#### Normailzed JSD
```
cd nyemtaay/tests/testdata/greekbears
nyemtaay -f bear_aln_sequences.fasta -m metadata_lat_long.csv -s LOC -calc norm-jsd -dm
```

#### Directionality 
```
nyemtaay -f bear_aln_sequences.fasta -m metadata_lat_long.csv -s LOC -calc directionality -dm
```

#### Partitioning based on normailzed JSD
```
nyemtaay -f bear_aln_sequences.fasta -m metadata_lat_long.csv -s LOC -calc norm-jsd-cluster -dm
```

#### Overlay of both Partitioning and Directionality  
```
nyemtaay -f bear_aln_sequences.fasta -m metadata_lat_long.csv -s LOC -calc directionality norm-jsd-cluster -dm
```


### Population Genetics

```
nyemtaay -f bear_aln_sequences.fasta -m metadata_lat_long.csv -s LOC -calc nei-fst
```
