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
    Name string
    // maps a chromosome name to a map of its sequences
    Chroms map[string](map[string]string)
}


func checkError(err error) {
    if err != nil {
        log.Fatal(err)
    }
}


func lines(str string) (numLines int, lines []string) {
    numLines = 0
    for i := 0; i < len(str); i++ {
        if str[i] == '\n' {
            numLines++
        }
    }
    lines = strings.Split(str, "\n")
    return numLines, lines
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


func parseGenome(genomeName string) RefGenome {
    var refGenome RefGenome
    refGenome.Name = genomeName
    chromFileInfos, err := ioutil.ReadDir(genomeName)
    checkError(err)
    // loop through every file in the reference directory
    for i := 0; i < len(chromFileInfos); i++ {
        chromFilename := strings.Join([]string{genomeName, chromFileInfos[i].Name()}, "/")
        fmt.Println("chromFilename:", chromFilename)
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align)
        if strings.HasSuffix(chromFilename, ".fa") {
            rawSeqBytes, err := ioutil.ReadFile(chromFilename)
            checkError(err)
            rawSeq := string(rawSeqBytes)
            numLines, seqLines := lines(rawSeq)
            for i := 0; i < numLines; i++ {
                seqLine := seqLines[i]
                if seqLine[0] == byte('>') {
                    fmt.Println(seqLine[1:])
                }
            }
        }
    }
    return refGenome
}


func main() {

    if len(os.Args) != 2 {
        fmt.Println(len(os.Args), "args supplied")
        fmt.Println("arg error - usage: ./minimize <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[1]
    refGenome := parseGenome(genomeName)
    fmt.Println(refGenome)

    rawRepeatsBytes, err := ioutil.ReadFile("dm3/dm3.fa.out")
    checkError(err)
    rawRepeats := string(rawRepeatsBytes)

    numLines, matchLines := lines(rawRepeats)
    fmt.Println("number of repeat file lines:", numLines)
    // drop header and empty end line
    matchLines = matchLines[3:len(matchLines)-1]
    fmt.Println("number of parsed matchLines:", len(matchLines))

    matches := parseMatchBares(matchLines)

    fmt.Println("number of parsed matches", len(matches))
}
