package leastsquares

import (
	"fmt"
	"log"
	"math"

	"gonum.org/v1/gonum/mat"
)

// Calc calculates the coefficients using the Least Squares method, using GoNum
func (rxn *Reaction) Calc() {
	cFormulas := make(map[string][]float64, len(rxn.UniqueElements))
	totalNumberOfCompounds := len(rxn.Reactants) + len(rxn.Products)
	matrixA := [][]float64{}
	matrixB := [][]float64{}

	for _, symbol := range rxn.UniqueElements {
		row := rxn.getUniqueMatrixRow(symbol)

		for i := range row {
			// If between len(rxn.Reactants) and totalNumberOfCompounds, switch the sign. (Like moving to other side of equation)
			if i > len(rxn.Reactants)-1 && i < totalNumberOfCompounds-1 && row[i] != 0 { // Don't flip 0 to -0
				row[i] = row[i] * -1
			}
		}
		cFormulas[symbol] = row
	}

	// Fill out Matrices A and B by splitting cFormulas into [everything but last column] and [last column], respectively
	for _, cFormulaRow := range cFormulas {
		tempMatrixARow := []float64{}
		tempMatrixBRow := []float64{}
		for c, cFormulaRowElement := range cFormulaRow {
			if c != totalNumberOfCompounds-1 {
				tempMatrixARow = append(tempMatrixARow, cFormulaRowElement)
			} else {
				tempMatrixBRow = append(tempMatrixBRow, cFormulaRowElement)
			}
		}
		matrixA = append(matrixA, tempMatrixARow)
		matrixB = append(matrixB, tempMatrixBRow)
	}

	/* printMatrix("A", matrixA)
	 * printMatrix("B", matrixB) */

	// Get my matrix data into a form that gonum can use.
	aData := []float64{}
	bData := []float64{}
	for _, row := range matrixA {
		for _, col := range row {
			aData = append(aData, col)
		}
	}
	for _, row := range matrixB {
		for _, col := range row {
			bData = append(bData, col)
		}
	}

	fmt.Println(matrixA)
	fmt.Println(matrixB)
	fmt.Printf("aData:%v\n", aData)
	fmt.Printf("bData:%v\n\n", bData)

	fmt.Printf("len(aData):%v\n", len(aData))
	fmt.Printf("len(bData):%v\n\n", len(bData))

	a := mat.NewDense(len(bData), len(bData), aData)
	b := mat.NewDense(len(bData), 1, bData)
	/* fmt.Printf("a:%v\n", a)
	 * fmt.Printf("b:%v\n\n", b) */

	// Compute the inverse of A.
	var aInv mat.Dense
	err := aInv.Inverse(a)
	var x mat.Dense
	if err != nil {
		log.Println("A is not invertible.  Using alternate method to solve.  See log for details.")
		/* log.Fatalf("A is not invertible: %v", err) */
		/* } */
		// Print the result using the formatter.
		/* fa := mat.Formatted(&aInv, mat.Prefix("       "), mat.Squeeze())
		 * fmt.Printf("aInv = %.2g\n\n", fa) */
		// The Inverse operation, however, should typically be avoided. If the
		// goal is to solve a linear system
		//  A * X = B,
		// then the inverse is not needed and computing the solution as
		// X = A^{-1} * B is slower and has worse stability properties than
		// solving the original problem. In this case, the SolveVec method of
		// VecDense (if B is a vector) or Solve method of Dense (if B is a
		// matrix) should be used instead of computing the Inverse of A.
		var tempX mat.Dense
		/* fmt.Println("[A] * [X] = [B] is a better, faster version of ")
		 * fmt.Println("[X] = [A]^{-1} * [B]") */
		err = tempX.Solve(a, b)
		if err != nil {
			log.Fatalf("no solution: %v", err)
		}
		x.Solve(a, b)

	} else {
		x.Mul(&aInv, b)
	}
	/* Print the result using the formatter. */
	/* fx := mat.Formatted(&x, mat.Prefix("    "), mat.Squeeze())
	 * fmt.Printf("X = %.1f = {coefficient 1, coefficient 2}\n", fx)
	 * fmt.Printf("\n%+v\n", x) */

	// Get all the coefficients but the last one
	coefficients := x.RawMatrix().Data

	lastCoefficient := math.Round(mat.Det(a))
	coefficients = append(coefficients, lastCoefficient)

	// Fix all signs to be positive
	for j := range coefficients {
		coefficients[j] = math.Abs(coefficients[j])
	}
	rxn.ComplexCoefficients = coefficients
	for i := range rxn.Compounds {
		rxn.Compounds[i].ComplexCoefficient = coefficients[i]
	}
}
