# # sabreur: fast, reliable and handy demultiplexing tool for fastx files.

## Download the data
We will use as example a demultiplex dataset freely available on figshare.

```bash
curl -L -o sabreur_data.tar.gz https://ndownloader.figshare.com/files/14481357
tar -xzvf sabreur_data.tar.gz
rm sabreur_data.tar.gz
cd sabreur_data
```

## Demultiplexing

sabreur is intended to be easy to use and handy. Therefore, as input it takes the barcode file and the input files without any other argument of options that need to be specified.

```
sabreur barcode.txt reads_R1.fq reads_R2.fq
```

### Compress the output
With sabreur you can easily compress the demultiplexed files even if the input is not. By default the demultiplexed files are compressed or not following the compression mode of the input.


#### Compress to gz
```
sabreur --format gz barcode.txt reads_R1.fq reads_R2.fq
```


#### Compress to xz
```
sabreur --format xz barcode.txt reads_R1.fq reads_R2.fq
```


#### Compress to bzip2
```
sabreur --format bz2 barcode.txt reads_R1.fq reads_R2.fq
```

It also possible to compress the demultiplexed files even if the input are not inputs.