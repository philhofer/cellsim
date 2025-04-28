package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"math/rand/v2"
	"os"
	"path/filepath"
	"slices"
)

func fatal(err error) {
	fmt.Fprintln(os.Stderr, err)
	os.Exit(1)
}

var (
	numCells    int
	numGenes    int
	umisPerCell int
	umisStddev  int
	numLowExp   int
	odir        string
)

func init() {
	flag.IntVar(&numCells, "cells", 20000, "number of cell barcodes")
	flag.IntVar(&numGenes, "genes", 400, "number of genes")
	flag.IntVar(&umisPerCell, "umis", 1000, "mean umis/cell")
	flag.IntVar(&umisStddev, "sigma", 250, "standard deviation in umis/cell")
	flag.IntVar(&numLowExp, "num-low-exp", 1, "number of low-expression genes")
	flag.StringVar(&odir, "o", ".", "output directory")
}

type gene struct {
	name   string
	weight int
}

type geneset struct {
	genes []gene
	cdf   []float64
}

// sample a single UMI and return the gene number
func (g *geneset) sample1() int {
	f := rand.Float64()
	n, _ := slices.BinarySearch(g.cdf, f)
	return n
}

func mkgeneset(n int, wgen func(int) int) geneset {
	lst := make([]gene, n)
	cdf := make([]float64, n)
	total := 0
	for i := range lst {
		w := wgen(i)
		total += w
		lst[i] = gene{
			name:   fmt.Sprintf("gene%d", i),
			weight: w,
		}
	}
	ft := float64(total)
	rsum := 0.0
	for i := range cdf {
		c := float64(lst[i].weight) / ft
		cdf[i] = c + rsum
		rsum = cdf[i]
	}
	return geneset{
		genes: lst,
		cdf:   cdf,
	}
}

func writeFeatures() {
	f, err := os.Create(filepath.Join(odir, "features.tsv.gz"))
	if err != nil {
		fatal(err)
	}
	defer f.Close()
	fgz := gzip.NewWriter(f)
	defer fgz.Close()
	for i := range numGenes {
		fmt.Fprintf(fgz, "gene%d\tgene%d\tGene Expression\n", i, i)
	}
}

var nts = [4]byte{'A', 'C', 'T', 'G'}

func barcode(size int) string {
	buf := make([]byte, size)
	for i := range buf {
		buf[i] = nts[rand.IntN(4)]
	}
	return string(buf)
}

func writeBarcodes() {
	f, err := os.Create(filepath.Join(odir, "barcodes.tsv.gz"))
	if err != nil {
		fatal(err)
	}
	defer f.Close()
	fgz := gzip.NewWriter(f)
	defer fgz.Close()
	for range numCells {
		fmt.Fprintf(fgz, "%s\n", barcode(16))
	}
}

func main() {
	flag.Parse()
	writeFeatures()
	writeBarcodes()
	matrix := make([]int, numGenes*numCells)

	const weight0 = 100
	gs := mkgeneset(numGenes, func(i int) int {
		if i < numLowExp {
			// relative rate of ~0.5
			return (numGenes * weight0) / (2 * umisPerCell)
		}
		return weight0
	})

	for cell := range numCells {
		size := int(rand.NormFloat64()*float64(umisStddev)) + umisPerCell
		size = max(size, 1)
		for range size {
			matrix[(cell*numGenes)+gs.sample1()]++
		}
	}
	nnz := 0
	for _, c := range matrix {
		if c != 0 {
			nnz++
		}
	}
	ofile, err := os.Create(filepath.Join(odir, "matrix.mtx.gz"))
	if err != nil {
		fatal(err)
	}
	defer ofile.Close()
	fgz := gzip.NewWriter(ofile)
	defer fgz.Close()
	fmt.Fprintln(fgz, "%%MatrixMarket matrix coordinate integer general")
	fmt.Fprintf(fgz, "%d %d %d\n", numGenes, numCells, nnz)
	for cell := range numCells {
		row := matrix[(cell * numGenes) : (cell+1)*numGenes]
		for gene, count := range row {
			if count == 0 {
				continue
			}
			fmt.Fprintf(fgz, "%d %d %d\n", gene+1, cell+1, count)
		}
	}
}
