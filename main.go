package main

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
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
		UnbalancedFormula: "Cr + O2 = Cr2O3",
	}
	reaction.Compounds = getCompoundDetails(reaction.UnbalancedFormula)
	reaction.Coefficients = calcCoefficients(reaction)
	fmt.Printf("Reaction Details:\n%s", prettify(reaction))

	// Have a list of all elements as variables
	// Piece together a formula using elements and subscripts as slices?
	/*     fmt.Println("Loading Element Symbols csv...")
	 *     symbCSV, err := ioutil.ReadFile("./Periodic_Table_of_Elements-Symbols.csv")
	 *     if err != nil {
	 *         log.Println(err)
	 *     }
	 *     symbString := string(symbCSV)
	 *     elementSymbolsCSV := csv.NewReader(strings.NewReader(symbString))
	 *
	 *     record, err := elementSymbolsCSV.Read()
	 *     if err != nil {
	 *         log.Fatal(err)
	 *     }
	 *
	 *     fmt.Println(record[0]) */

}

// Reaction is a chemical reaciton to be analyzed
type Reaction struct {
	UnbalancedFormula string     `json:"unbalanced_formula"`
	BalancedFormula   string     `json:"balanced_formula"`
	Coefficients      []float64  `json:"coefficients"`
	Compounds         []Compound `json:"compounds"`
	Reactants         []Compound `json:"reactants"`
	Products          []Compound `json:"products"`
}

// Compound is a set of elements joined together in a chemical reaction
type Compound struct {
	Formula     string    `json:"formula"`
	Coefficient float64   `json:"coefficient"`
	Elements    []Element `json:"elements"`
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

func getCompoundDetails(formula string) []Compound {
	r := regexp.MustCompile(`([a-zA-z]{1,2}\d?)*`)
	cs := []string{}
	compounds := []Compound{}
	c := r.FindAllString(formula, -1)
	for e := range c {
		if c[e] != "" {
			cs = append(cs, c[e])
		}
	}
	for _, compound := range cs {
		tempCompound := Compound{
			Formula: compound,
		}
		tempCompound.Elements = getElementsInCompound(tempCompound)
		compounds = append(compounds, tempCompound)
	}
	return compounds
}

func getElementsInCompound(compound Compound) []Element {
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

	}
	return elements
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

func calcCoefficients(rxn Reaction) []float64 {
	var placeholderFloat64 []float64 = []float64{100.}
	return placeholderFloat64
}
