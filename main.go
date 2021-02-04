package main

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"
)

func init() {
	// Create log file
	// If the file doesn't exist, create it or append to the file
	file, err := os.OpenFile("logs.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	log.SetOutput(file)
}

// Reaction is a chemical reaciton to be analyzed
type Reaction struct {
	UnbalancedFormula      string     `json:"unbalanced_formula"`
	BalancedFormula        string     `json:"balanced_formula"`
	ComplexCoefficients    []float64  `json:"complex_coefficients"`
	SimplifiedCoefficients []float64  `json:"simplified_coefficients"`
	Compounds              []Compound `json:"compounds"`
	Reactants              []string   `json:"reactants"`
	Products               []string   `json:"products"`
	UniqueElements         []string   `json:"unique_elements"`
}

// Compound is a set of elements joined together in a chemical reaction
type Compound struct {
	Formula               string    `json:"formula"`
	ComplexCoefficient    float64   `json:"complex_coefficient"`
	SimplifiedCoefficient float64   `json:"simplified_coefficient"`
	Elements              []Element `json:"elements"`
	Index                 []int     `json:"index"`
	Reactant              bool      `json:"reactant"`
	Product               bool      `json:"product"`
}

// Element is an individual part of a compound
type Element struct {
	Symbol    string  `json:"symbol"`
	Subscript float64 `json:"subscript"`
}

func prettify(i interface{}) string {
	s, _ := json.MarshalIndent(i, "", "  ")
	return string(s)
}

// Method to get details of all the compounds in the reaction.
func (rxn *Reaction) getCompoundDetails() {
	formula := rxn.UnbalancedFormula

	r := regexp.MustCompile(`([a-zA-z]{1,2}\d?)*`)
	r2 := regexp.MustCompile(`[=]|[<]?[-]{1,3}[>]?`) // Match either = or (<)--> for yields symbol
	cs := []string{}
	csIndex := [][]int{}

	// Find location of yields symbol.  Returns an array of two values, the first is the starting index, and the second is the non-inclusive ending index.
	yieldsIndex := r2.FindStringIndex(formula)
	/* fmt.Printf("The index of the yields symbol is at %d\n", yieldsIndex) */

	c := r.FindAllString(formula, -1)
	cIndex := r.FindAllStringIndex(formula, -1)
	for e := range c {
		if c[e] != "" {
			cs = append(cs, c[e])
			csIndex = append(csIndex, cIndex[e])
		}
	}
	for i, compound := range cs {
		tempCompound := Compound{
			Formula: compound,
			Index:   csIndex[i],
		}
		// If the compound's index comes before the yields symbol, it is a reactant, otherwise, it is a product.
		if tempCompound.Index[0] < yieldsIndex[0] {
			tempCompound.Reactant = true
			rxn.Reactants = append(rxn.Reactants, tempCompound.Formula)
		} else {
			tempCompound.Product = true
			rxn.Products = append(rxn.Products, tempCompound.Formula)
		}
		tempCompound.Elements = rxn.getElementsInCompound(tempCompound)
		rxn.Compounds = append(rxn.Compounds, tempCompound)
	}
}

// Method to get the Elements in a compound of a reaction, and update the reaction.  It returns a slice of Elements
func (rxn *Reaction) getElementsInCompound(compound Compound) []Element {
	r := regexp.MustCompile(`[A-Z]{1}([a-z]?)`)
	elements := []Element{}
	es := []string{}

	e := r.FindAllString(compound.Formula, -1)
	// e is a []string, so we loop through it to pull out each individual element, discarding the blanks
	for elem := range e {
		if e[elem] != "" {
			/* fmt.Printf("Found element in %s: %s\n", compound.Formula, e[elem]) */
			es = append(es, e[elem])
		}
	}
	// Now we have one slice containing all of the elements as strings (es)

	// Get a matching slice of the subscripts in the compound
	subscripts := getSubscriptsInCompound(compound.Formula)
	for e, element := range es {
		// Convert to an Element struct with only the formula filled out, and add to the slice
		tempElement := Element{
			Symbol:    element,
			Subscript: subscripts[e],
		}
		elements = append(elements, tempElement)

		// Check to see if this element is already added to the unique elements array
		isInUnique := containsString(rxn.UniqueElements, tempElement.Symbol)
		if isInUnique == false {
			rxn.UniqueElements = append(rxn.UniqueElements, tempElement.Symbol)
		}
	}
	return elements
}

func containsString(s []string, e string) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func getSubscriptsInCompound(formula string) []float64 {
	r := regexp.MustCompile(`([A-Z][a-z]?)(\d*)`)
	s := r.FindAllStringSubmatch(formula, -1)
	subscripts := []float64{}
	for _, sub := range s {
		if sub[2] == "" {
			sub[2] = "1"
		}
		subf64, _ := strconv.ParseFloat(sub[2], 64)
		subscripts = append(subscripts, subf64)
	}

	return subscripts
}

func (rxn *Reaction) calcComplexCoefficients() {
	cFormulas := make(map[string][]float64, len(rxn.UniqueElements))
	totalNumberOfCompounds := len(rxn.Reactants) + len(rxn.Products)
	matrixA := [][]float64{}
	matrixB := [][]float64{}

	for _, s := range rxn.UniqueElements {
		/* row := []float64{1, 1, 2, 3} */
		row := rxn.getUniqueMatrixRow(s)

		for i := range row {
			// If between len(rxn.Reactants) and totalNumberOfCompounds, switch the sign. (Like moving to other side of equation)
			if i > len(rxn.Reactants)-1 && i < totalNumberOfCompounds-1 {
				// Don't flip 0 to -0
				if row[i] != 0 {
					row[i] = row[i] * -1
				}
			}
		}
		cFormulas[s] = row
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

	/* fmt.Println(matrixA)
	 * fmt.Println(matrixB)
	 * fmt.Printf("aData:%v\n", aData)
	 * fmt.Printf("bData:%v\n\n", bData) */

	a := mat.NewDense(len(matrixA[0]), len(matrixA[0]), aData)
	b := mat.NewDense(len(bData), 1, bData)
	// Initialize a matrix A.
	/*     a := mat.NewDense(3, 3, []float64{
	 *         1, 0, 0,
	 *         1, 0, -3,
	 *         0, 1, -2,
	 *     })
	 *
	 *     b := mat.NewDense(3, 1, []float64{
	 *         1,
	 *         0,
	 *         0,
	 *     }) */
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
	// Print the result using the formatter.
	/* fx := mat.Formatted(&x, mat.Prefix("    "), mat.Squeeze())
	 * fmt.Printf("X = %.1f = {coefficient 1, coefficient 2}\n", fx)
	 * fmt.Printf("\n%+v\n", x) */

	// Get all the coefficients but the last one
	coefficients := x.RawMatrix().Data

	lastCoefficient := math.Round(mat.Det(a))
	coefficients = append(coefficients, lastCoefficient)
	fmt.Printf("Unsimplified Coefficients: %v\n", coefficients)
	rxn.ComplexCoefficients = coefficients
	for i := range rxn.Compounds {
		rxn.Compounds[i].ComplexCoefficient = coefficients[i]
	}
}

func (rxn *Reaction) calcSimplifiedCoefficients() {
	vals := rxn.ComplexCoefficients
	tableOfVals := makeTable()
	fractionVals := [][]float64{}

	for _, val := range vals {
		for _, tableVal := range tableOfVals {
			if val == tableVal[0] {
				fractionVals = append(fractionVals, tableVal)
				break
			}
		}
	}

	/* fmt.Printf("fractionVals: %v\n", fractionVals) */
	//fmt.Printf("f: %v\n\n", f)

	numerators := []float64{}
	denominators := []float64{}

	for i := range fractionVals {
		numerators = append(numerators, fractionVals[i][1])
		denominators = append(denominators, fractionVals[i][2])
	}

	for d, denom := range denominators {
		if denom != 1 {

			for n, numer := range numerators {
				numerators[n] = (numer * denom)
			}
			numerators[d] /= denom
			denominators[d] = 1
			//break
		}

		values := [][]float64{}
		for j := range numerators {
			values = append(values, []float64{numerators[j], denominators[j]})
		}
		fmt.Println(values)
	}

	gcf := numerators[0]
	for _, a := range numerators {
		gcf = calcGCD(gcf, a)
	}

	for b := range numerators {
		numerators[b] = numerators[b] / gcf
	}

	rxn.SimplifiedCoefficients = numerators
}

func makeTable() [][]float64 {
	values := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25}
	decimals := [][]float64{}
	for _, n := range values {
		for _, d := range values {
			decimals = append(decimals, []float64{n / d, n, d})
		}
	}
	return decimals
}

func calcGCD(a, b float64) float64 {
	if b == 0 {
		return a
	}
	return calcGCD(b, math.Mod(a, b))
}

func getSymbolsFromCSV() []string {
	// Have a list of all elements as variables
	// Piece together a formula using elements and subscripts as slices?
	fmt.Println("Loading Element Symbols csv...")
	symbCSV, err := ioutil.ReadFile("./Periodic_Table_of_Elements-Symbols.csv")
	if err != nil {
		log.Println(err)
	}
	symbString := string(symbCSV)
	elementSymbolsCSV := csv.NewReader(strings.NewReader(symbString))

	record, err := elementSymbolsCSV.Read()
	if err != nil {
		log.Fatal(err)
	}

	return record
}

func printMatrix(name string, matrix [][]float64) {
	fmt.Printf("\n[%s]:\n", name)
	for _, row := range matrix {
		fmt.Println(row)
	}
	fmt.Println("")
}

// Generates a unique row of a matrix by comparing each unique element in the formula with how many times it shows up in each compound.
func (rxn *Reaction) getUniqueMatrixRow(symbol string) []float64 {
	matrixRow := []float64{}
	for _, compound := range rxn.Compounds {
		hasSymbol := false
		var subscript float64
		for _, element := range compound.Elements {
			if element.Symbol == symbol {
				hasSymbol = true
				subscript = element.Subscript
			}
		}
		if hasSymbol == true {
			matrixRow = append(matrixRow, subscript)
		} else {
			matrixRow = append(matrixRow, float64(0))
		}
	}
	return matrixRow
}

func (rxn *Reaction) createBalancedRxn() {
	var balancedRxn strings.Builder
	// Append reactants with coefficients
	for i := range rxn.Reactants {
		if i != 0 {
			balancedRxn.WriteString(" + ")
		}
		coefficientFloat := strconv.FormatFloat(rxn.SimplifiedCoefficients[i], 'f', -1, 64)
		balancedRxn.WriteString(coefficientFloat)
		balancedRxn.WriteString(rxn.Reactants[i])
	}

	// Add yields symbol
	balancedRxn.WriteString(" = ")

	// Append products with coefficients
	for i := range rxn.Products {
		if i != 0 {
			balancedRxn.WriteString(" + ")
		}
		coefficientFloat := strconv.FormatFloat(rxn.SimplifiedCoefficients[i+len(rxn.Reactants)], 'f', -1, 64)
		balancedRxn.WriteString(coefficientFloat)
		balancedRxn.WriteString(rxn.Products[i])
	}
	rxn.BalancedFormula = balancedRxn.String()
}

func main() {
	reaction := Reaction{
		/* UnbalancedFormula: "Cr + O2 = Cr2O3", */
		UnbalancedFormula: "MgO + Fe = Mg + Fe2O3",
		/* UnbalancedFormula: "NaHCO3 = Na2CO3 + CO2 + H2O", */
	}
	reaction.getCompoundDetails()
	reaction.calcComplexCoefficients()
	reaction.calcSimplifiedCoefficients()
	reaction.createBalancedRxn()
	fmt.Printf("%s\n", prettify(reaction))
}
