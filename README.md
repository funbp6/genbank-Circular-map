# Circular map
circular map plot generator, with genbank file input  
## Usage  
This plotting program can run with default setting  
input your genbank file, and svg circular map will be generated.   
```
./COGplot.py GENBANK.gb  
```
If you only want to draw with partial contigs 
```
./COGplot.py GENBANK.gb -p LIST_OF_CONTIGS_NAME  
```
other detail setting please check by `./COGplot.py -h` command  
## library required  
**Python3**  
* Biopython   
* colour  
* svgwrite  
* argparse  
## Example

![image](https://github.com/funbp6/genbank-Circular-map/blob/master/circular_map1.png)



