package main

import (
	"encoding/csv"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"regexp"
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

	fmt.Println(record[0])

	formula := "Cr + O2 --> Cr2O3"
	fmt.Printf("Formula: %s\n", formula)
	compounds := getCompounds(formula)
	fmt.Printf("Compounds: %v\n", compounds)
	for i := range compounds {
		elements := getElementsInCompound(compounds[i])
		fmt.Printf("\nElements in %v: %v\n", compounds[i], elements)
		subscripts := getSubscriptsInCompound(compounds[i])
		fmt.Printf("Subscript(s) for %v: %v\n", compounds[i], subscripts)
	}
}

func getCompounds(formula string) []string {
	r := regexp.MustCompile(`([a-zA-z]{1,2}\d?)*`)
	compounds := []string{}
	c := r.FindAllString(formula, -1)
	for e := range c {
		if c[e] != "" {
			compounds = append(compounds, c[e])
		}
	}
	return compounds
}

func getElementsInCompound(compound string) []string {
	r := regexp.MustCompile(`[A-Z]{1}([a-z]?)`)
	elements := []string{}
	es := r.FindAllString(compound, -1)
	for e := range es {
		if es[e] != "" {
			elements = append(elements, es[e])
		}
	}
	return elements
}

func getSubscriptsInCompound(compound string) [][]string {
	r := regexp.MustCompile(`([A-Z][a-z]?)(\d*)`)
	subscripts := r.FindAllStringSubmatch(compound, -1)
	for _, subscript := range subscripts {
		if subscript[2] == "" {
			subscript[2] = "1"
		}
	}
	return subscripts
}
