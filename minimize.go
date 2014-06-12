// A barebones (at the moment) Go script for parsing and minimizing repeats

package main

import ("fmt"
        "log"
        "io/ioutil"
        "strings"
)

func main() {

    var err error
    rawRepeatsBytes, err := ioutil.ReadFile("dm3/dm3.fa.out")
    rawRepeats := string(rawRepeatsBytes)
    if err != nil {
        log.Fatal(err)
    }

    numLines := 0
    for i := 0; i < len(rawRepeats); i++ {
        if rawRepeats[i] == '\n' {
            numLines++
        }
    }
    fmt.Println("number of repeat file lines:", numLines)

    repeats := strings.Split(rawRepeats, "\n")
    // drop header
    repeats = repeats[3:len(repeats)-1]
    fmt.Println("number of parsed repeats:", len(repeats))
    fmt.Println("first", repeats[0])
    fmt.Println("second", repeats[1])
    fmt.Println("third", repeats[2])
    fmt.Println("last", repeats[len(repeats)-1])
}
