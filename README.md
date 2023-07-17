## Usages

```
./prepare_isoform_data.py -g genome.fa -i circRNAs.isoforms.tsv -o out
```

where the `circRNAs.isoforms.tsv` is consist of the isoform_IDs,
e.g.
```
10|100036449|100036604|+
10|100042282,100048758,100054347,100057013,100063614,100065188|100042573,100048876,100054446,100057152,100063725,100065309|-
```

### Dependency
- Python 3
- bedtools
