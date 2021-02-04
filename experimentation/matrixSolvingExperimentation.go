package main

import (
	"fmt"
	"log"

	"gonum.org/v1/gonum/mat"
)

type matrix [][]float64

func (m matrix) print() {
	for _, r := range m {
		fmt.Println(r)
	}
	fmt.Println("")
}

func printSymMatrix(x mat.SymDense) {
	fx := mat.Formatted(&x, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("x = %v\n\n", fx)
}

func printDenseMatrix(x mat.Dense) {
	fx := mat.Formatted(&x, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("x = %v\n\n", fx)
}

func calcRREF(m matrix) {
	lead := 0
	rowCount := len(m)
	columnCount := len(m[0])
	for r := 0; r < rowCount; r++ {
		if lead >= columnCount {
			return
		}
		i := r
		for m[i][lead] == 0 {
			i++
			if rowCount == i {
				i = r
				lead++
				if columnCount == lead {
					return
				}
			}
		}
		m[i], m[r] = m[r], m[i]
		f := 1 / m[r][lead]
		for j := range m[r] {
			m[r][j] *= f
		}
		for i = 0; i < rowCount; i++ {
			if i != r {
				f = m[i][lead]
				for j, e := range m[r] {
					m[i][j] -= e * f
				}
			}
		}
		lead++
	}

	m.print()
}

func format(vals []float64, prec int, eps float64) []string {
	s := make([]string, len(vals))
	for i, v := range vals {
		if v < eps {
			s[i] = fmt.Sprintf("<%.*g", prec, eps)
			continue
		}
		s[i] = fmt.Sprintf("%.*g", prec, v)
	}
	return s
}

func calcSVD(a, b mat.Matrix) {
	// Perform an SVD retaining all singular vectors.
	var svd mat.SVD
	ok := svd.Factorize(a, mat.SVDFull)
	if !ok {
		log.Fatal("failed to factorize A")
	}

	// Determine the rank of the A matrix with a near zero condition threshold.
	const rcond = 1e-15
	rank := svd.Rank(rcond)
	if rank == 0 {
		log.Fatal("zero rank system")
	}

	// Find a least-squares solution using the determined parts of the system.
	var x mat.Dense
	svd.SolveTo(&x, b, rank)

	fmt.Printf("singular values = %v\nrank = %d\nx = %.15f",
		format(svd.Values(nil), 4, rcond), rank, mat.Formatted(&x, mat.Prefix("    ")))
}

func main() {
	/* "Cr + O2 = Cr2O3" */
	m := matrix{
		{1, 0, 2},
		{0, 2, 3},
	}

	/* "Fe + AgNO3 = Fe(NO3)2 + Ag", */
	/* m := matrix{
	 *     {1, 0, -1, 0},
	 *     {0, 1, 0, 1},
	 *     {0, 1, 2, 0},
	 *     {0, 3, 6, 0},
	 * } */

	/* "Na2CO3 + CO2 + H2O = NaHCO3", */
	/* m := matrix{
	 *     {2, 0, 0, 2},
	 *     {1, 1, 0, 3},
	 *     {3, 2, 1, 3},
	 *     {0, 0, 2, 0},
	 * } */

	/* "Na2CO3 + CO2 + H2O = NaHCO3", */
	/*     a := mat.NewDense(4, 3, []float64{
	 *         2, 0, 0,
	 *         1, 1, 0,
	 *         3, 2, 1,
	 *         0, 0, 2,
	 *     })
	 *
	 *     b := mat.NewDense(4, 1, []float64{
	 *         2,
	 *         3,
	 *         3,
	 *         0,
	 *     })
	 *     calcSVD(a, b) */

	calcRREF(m)

}
