package main

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"
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

func main() {
	reaction := Reaction{
		/* UnbalancedFormula: "Cr + O2 = Cr2O3", */
		/* UnbalancedFormula: "MgO + Fe = Fe2O3 + Mg", */
		UnbalancedFormula: "NaHCO3 = Na2CO3 + CO2 + H2O",
	}
	reaction.getCompoundDetails()
	reaction.calcCoefficients()
	fmt.Printf("%s\n", prettify(reaction))

}

// Reaction is a chemical reaciton to be analyzed
type Reaction struct {
	UnbalancedFormula string     `json:"unbalanced_formula"`
	BalancedFormula   string     `json:"balanced_formula"`
	Coefficients      []float64  `json:"coefficients"`
	Compounds         []Compound `json:"compounds"`
	Reactants         []string   `json:"reactants"`
	Products          []string   `json:"products"`
	UniqueElements    []string   `json:"unique_elements"`
}

// Compound is a set of elements joined together in a chemical reaction
type Compound struct {
	Formula     string    `json:"formula"`
	Coefficient float64   `json:"coefficient"`
	Elements    []Element `json:"elements"`
	Index       []int     `json:"index"`
	Reactant    bool      `json:"reactant"`
	Product     bool      `json:"product"`
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

func (rxn *Reaction) calcCoefficients() {
	cFormulas := make(map[string][]float64, len(rxn.UniqueElements))
	totalNumberOfCompounds := len(rxn.Reactants) + len(rxn.Products)
	matrixA := [][]float64{}
	matrixB := [][]float64{}

	for _, s := range rxn.UniqueElements {
		row := []float64{1, 1, 2, 3}

		for i := range row {
			// If between len(rxn.Reactants) and totalNumberOfCompounds, switch the sign. (Like moving to other side of equation)
			if i > len(rxn.Reactants)-1 && i < totalNumberOfCompounds-1 {
				row[i] = row[i] * -1
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

    printMatrix("A", matrixA)
    printMatrix("B", matrixB)
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
