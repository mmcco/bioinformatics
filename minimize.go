// A barebones (at the moment) Go script for parsing and minimizing repeats

package main


import ("fmt"
        "log"
        "os"
        "io/ioutil"
        "strings"
        "strconv"
)

var err error

// uses 64-bit values, which is probably unnecessary for some fields
// this is done for simplicity (and strconv compatibility)
// it should probably be changed eventually
type Match struct {
    SW_Score int64
    PercDiv float64
    PercDel float64
    PercIns float64
    SeqName string
    SeqStart int64
    SeqEnd int64
    Left int64
    IsComplement bool
    RepeatType string
    RepeatClass string
    HangingBases int64
    RepeatStart int64
    RepeatEnd int64
    ID int64
}


// a minimal Match implementation for minimizing
type MatchBare struct {
    SeqName string
    SeqStart int64
    SeqEnd int64
    IsComplement bool
    RepeatType string
    RepeatClass string
    ID int64
}


type RefGenome struct {
    Chroms []string
}


func checkError(err error) {
    if err != nil {
        log.Fatal(err)
    }
}

func parseMatches(matchLines []string) []Match {
    var matches []Match
    for i := 0; i < len(matchLines); i++ {
        rawVals := strings.Fields(matchLines[i])
        var match Match
        match.SW_Score, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.PercDiv, err = strconv.ParseFloat(rawVals[1], 64)
        checkError(err)
        match.PercDel, err = strconv.ParseFloat(rawVals[2], 64)
        checkError(err)
        match.PercIns, err = strconv.ParseFloat(rawVals[3], 64)
        checkError(err)
        match.SeqName = rawVals[4]
        match.SeqStart, err = strconv.ParseInt(rawVals[5], 10, 64)
        checkError(err)
        match.SeqEnd, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.Left, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.IsComplement = rawVals[8] == "C"
        match.RepeatType = rawVals[9]
        match.RepeatClass = rawVals[10]
        match.HangingBases, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.RepeatStart, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.RepeatEnd, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.ID, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)

        matches = append(matches, match)
    }
    return matches
}


func parseMatchBares(matchLines []string) []MatchBare {
    var matches []MatchBare
    for i := 0; i < len(matchLines); i++ {
        rawVals := strings.Fields(matchLines[i])
        var match MatchBare
        match.SeqName = rawVals[4]
        match.SeqStart, err = strconv.ParseInt(rawVals[5], 10, 64)
        checkError(err)
        match.SeqEnd, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.IsComplement = rawVals[8] == "C"
        match.RepeatType = rawVals[9]
        match.RepeatClass = rawVals[10]
        match.ID, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)

        matches = append(matches, match)
    }
    return matches
}


func main() {

    if len(os.Args) != 2 {
        fmt.Println(len(os.Args), "args supplied")
        fmt.Println("arg error - usage: ./minimize <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[1]
    fmt.Println(genomeName)

    rawRepeatsBytes, err := ioutil.ReadFile("dm3/dm3.fa.out")
    rawRepeats := string(rawRepeatsBytes)
    checkError(err)

    numLines := 0
    for i := 0; i < len(rawRepeats); i++ {
        if rawRepeats[i] == '\n' {
            numLines++
        }
    }
    fmt.Println("number of repeat file lines:", numLines)

    matchLines := strings.Split(rawRepeats, "\n")
    // drop header and empty end line
    matchLines = matchLines[3:len(matchLines)-1]
    fmt.Println("number of parsed matchLines:", len(matchLines))

    matches := parseMatchBares(matchLines)

    fmt.Println("number of parsed matches", len(matches))
}
