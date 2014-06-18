/*
    WARNING!!! This program is currently under development and may be buggy or broken.

    A barebones (at the moment) Go script for parsing and minimizing repeats

    The sole command line argument is the name of the reference genome (e.g. "dm3").

    This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
        * a RepeatMasker library containing:
            - the match library (e.g. "dm3.fa.out")
            - the alignment information (e.g. "dm3.fa.align")
        * one or more reference genome files in FASTA format with the ".fa" filetype
*/


package main


import ("fmt"
        "log"
        "os"
        "io/ioutil"
        "strings"
        "strconv"
        "reflect"
        //"math"
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
    SeqRemains int64
    IsRevComp bool
    RepeatName string
    RepeatClass string
    RepeatStart int64
    RepeatEnd int64
    RepeatRemains int64
    RepeatID int64
}


// a minimal Match implementation for minimizing
type MatchBare struct {
    SeqName string
    SeqStart int64
    SeqEnd int64
    IsRevComp bool
    RepeatName string
    RepeatClass string
    RepeatID int64
}


type RefGenome struct {
    Name string
    // maps a chromosome name to a map of its sequences
    Chroms map[string](map[string]string)
}


type Kmer struct {
    Kmer string
    MinOffset int
    IsRevComp bool
    LCA *Repeat
}


type Repeat struct {
    ID int64
    Class []string
    FullName string
}


type ClassNode struct {
    Name string
    Class []string
    Parent *ClassNode
    Children []*ClassNode
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
    // drop the trailing newline's line if it's there
    if strings.TrimSpace(lines[len(lines)-1]) == "" {
        lines = lines[:len(lines)-1]
    }
    return numLines, lines
}


func parseMatches(matchLines []string) []Match {
    var matches []Match
    for i := 0; i < len(matchLines); i++ {
        rawVals := strings.Fields(matchLines[i])
        var match Match
        match.IsRevComp = rawVals[8] == "C"

        // remove enclosing parentheses
        // !!! in the future, checkes to ensure that the parentheses exist should be added
        // !!! it would also be sensible to check that rawVals[8] is either "C" or "+"
        rawVals[7] = rawVals[7][1:len(rawVals[7])-1]
        if match.IsRevComp {
            rawVals[11] = rawVals[11][1:len(rawVals[11])-1]
        } else {
            rawVals[13] = rawVals[13][1:len(rawVals[13])-1]
        }

        match.SW_Score, err = strconv.ParseInt(rawVals[0], 10, 64)
        checkError(err)
        match.PercDiv, err = strconv.ParseFloat(rawVals[1], 64)
        checkError(err)
        match.PercDel, err = strconv.ParseFloat(rawVals[2], 64)
        checkError(err)
        match.PercIns, err = strconv.ParseFloat(rawVals[3], 64)
        checkError(err)
        match.SeqName = strings.TrimSpace(rawVals[4])
        match.SeqStart, err = strconv.ParseInt(rawVals[5], 10, 64)
        checkError(err)
        match.SeqEnd, err = strconv.ParseInt(rawVals[6], 10, 64)
        checkError(err)
        match.SeqRemains, err = strconv.ParseInt(rawVals[7], 10, 64)
        checkError(err)
        // match.IsComplement, rawVals[8], moved above
        match.RepeatName = strings.TrimSpace(rawVals[9])
        match.RepeatClass = strings.TrimSpace(rawVals[10])
        match.RepeatStart, err = strconv.ParseInt(rawVals[11], 10, 64)
        checkError(err)
        match.RepeatEnd, err = strconv.ParseInt(rawVals[12], 10, 64)
        checkError(err)
        match.RepeatRemains, err = strconv.ParseInt(rawVals[13], 10, 64)
        checkError(err)
        match.RepeatID, err = strconv.ParseInt(rawVals[14], 10, 64)
        checkError(err)

        // necessary swaps to convert reverse complement repeat indexes to positive-strand indexes
        if match.IsRevComp {
            match.RepeatStart, match.RepeatEnd, match.RepeatRemains =
                match.RepeatRemains, match.RepeatStart, match.RepeatRemains + (match.RepeatEnd - match.RepeatStart)
        }

        matches = append(matches, match)
    }
    return matches
}


func parseMatchBares(matchLines []string) []MatchBare {
    var matches []MatchBare
    for i := 0; i < len(matchLines); i++ {
        rawVals := strings.Fields(matchLines[i])
        var match MatchBare
        match.SeqName = strings.TrimSpace(rawVals[4])
        match.SeqStart, err = strconv.ParseInt(rawVals[5], 10, 64)
        checkError(err)
        match.SeqEnd, err = strconv.ParseInt(rawVals[6], 10, 64)
        checkError(err)
        match.IsRevComp = rawVals[8] == "C"
        match.RepeatName = strings.TrimSpace(rawVals[9])
        match.RepeatClass = strings.TrimSpace(rawVals[10])
        match.RepeatID, err = strconv.ParseInt(rawVals[14], 10, 64)
        checkError(err)

        matches = append(matches, match)
    }
    return matches
}


func parseGenome(genomeName string) RefGenome {
    refGenome := RefGenome{genomeName, make(map[string](map[string]string))}
    chromFileInfos, err := ioutil.ReadDir(genomeName)
    checkError(err)
    // used below to store the two keys for RefGenome.Chroms
    for i := 0; i < len(chromFileInfos); i++ {
        chromFilename := strings.Join([]string{genomeName, chromFileInfos[i].Name()}, "/")
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align)
        if strings.HasSuffix(chromFilename, ".fa") {
            rawSeqBytes, err := ioutil.ReadFile(chromFilename)
            checkError(err)
            rawSeq := string(rawSeqBytes)
            numLines, seqLines := lines(rawSeq)

            // ultimately contains this file's sequences
            seqMap := make(map[string]string)

            // populate thisSeq with the first seq's name
            thisSeq := append([]string(nil), strings.TrimSpace(seqLines[0])[1:])
            for i := 1; i < numLines; i++ {
                seqLine := strings.TrimSpace(seqLines[i])
                if seqLine[0] == byte('>') {
                    // we now have a full seq, and can write it to the map
                    seqMap[thisSeq[0]] = strings.Join(thisSeq[1:], "")
                    thisSeq = []string{seqLine[1:]}
                } else {
                    thisSeq = append(thisSeq, seqLine)
                }
            }
            // add the remaining sequence
            seqMap[thisSeq[0]] = strings.Join(thisSeq[1:], "")
            // finally, we insert this map into the full map
            chromName := chromFilename[len(genomeName)+1:len(chromFilename)-3]
            // must initialize the inner map
            refGenome.Chroms[chromName] = make(map[string]string)
            refGenome.Chroms[chromName] = seqMap
        }
    }
    return refGenome
}


func getMinimizer(kmer string, m int) int {
    currMin := 0
    for i := 0; i < len(kmer) - m + 1; i++ {
        if kmer[i:i+m] < kmer[currMin:currMin+m] {
            currMin = i
        }
    }
    return currMin
}


func reverseComplement(seq string) string {
    var revCompSeq []byte
    for i := 0; i < len(seq); i++ {
            switch seq[len(seq) - i - 1] {
            case 'a': revCompSeq = append(revCompSeq, 't')
            case 'A': revCompSeq = append(revCompSeq, 'T')
            case 't': revCompSeq = append(revCompSeq, 'a')
            case 'T': revCompSeq = append(revCompSeq, 'A')
            case 'c': revCompSeq = append(revCompSeq, 'g')
            case 'C': revCompSeq = append(revCompSeq, 'G')
            case 'g': revCompSeq = append(revCompSeq, 'c')
            case 'G': revCompSeq = append(revCompSeq, 'C')
            }
    }
    return string(revCompSeq)
}


func max(a int64, b int64) int64 {
    if a > b {
        return a
    } else {
        return b
    }
}


func min(a, b int) int {
    if a < b {
        return a
    } else {
        return b
    }
}


func (repeat Repeat) Print() {
    for j := range repeat.Class {
        for j_ := 0; j_ < j; j_++ {
            fmt.Printf("\t")
        }
        fmt.Printf("%s\n", repeat.Class[j])
    }
    fmt.Println()
}


func (classNode *ClassNode) PrintTree() {
    classNode.printTreeRec(0)
}

func (classNode *ClassNode) printTreeRec(indent int) {
    for i := 0; i < indent; i++ {
        fmt.Printf("\t")
    }
    if len(classNode.Class) == 0 {
        fmt.Println("root")
    } else {
        fmt.Println(classNode.Name)
        //fmt.Printf("%s\n", classNode.Class[indent])
    }
    for i := range classNode.Children {
        classNode.Children[i].printTreeRec(indent+1)
    }
}


func classSliceContains(a_s []*ClassNode, b *ClassNode) bool {
    for i := range a_s {
        if a_s[i] == b {
            return true
        }
    }
    return false
}

// returns a pointer to the root of the tree a map of ClassNode names to ClassNode pointers
func GetClassTree(repeats []Repeat) (*ClassNode, map[string](*ClassNode)) {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    classNodes := make(map[string](*ClassNode))
    // would be prettier if expanded
    classRoot := &ClassNode{}
    classNodes["root"] = classRoot
    classRoot.Name = "root"
    // all but Name is left nil
    for i := 1; i < len(repeats); i++ {
        // ignore the null indices
        if repeats[i].ID != 0 {
            // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
            if repeats[i].FullName == "root" {
                fmt.Println(repeats[i])
            }
            for j := 1; j <= len(repeats[i].Class); j++ {
                thisClass := repeats[i].Class[:j]
                thisClassName := strings.Join(thisClass, "/")
                var keyExists bool
                _, keyExists = classNodes[thisClassName]
                if !keyExists {
                    thisClassNode := ClassNode{}
                    classNodes[thisClassName] = &thisClassNode
                    thisClassNode.Name = thisClassName
                    classNodes[thisClassName].Class = thisClass
                    // first case handles primary classes, as classRoot is implicit and not listed in thisClass
                    if j == 1 {
                        classNodes[thisClassName].Parent = classRoot
                    } else {
                        classNodes[thisClassName].Parent = classNodes[strings.Join(thisClass[:len(thisClass)-1], "/")]
                    }
                    parent := classNodes[thisClassName].Parent
                    if parent.Children == nil {
                        parent.Children = []*ClassNode{}
                    }
                    parent.Children = append(parent.Children, classNodes[thisClassName])
                }
            }
        }
    }
    return classRoot, classNodes
}


func GetLCA(classNodes map[string](*ClassNode), cnA, cnB *ClassNode) *ClassNode {
    for i := min(len(cnA.Class), len(cnB.Class)); i > 0; i-- {
        if reflect.DeepEqual(cnA.Class[:i], cnB.Class[:i]) {
            // we walk back to the LCA's pointer through the parent fields
            lca := cnA
            for j := 0; j < len(cnA.Class) - i; i++ {
                lca = lca.Parent
            }
            return lca
        }
    }
    return classNodes["root"]
}


func main() {

    if len(os.Args) != 2 {
        fmt.Println(len(os.Args), "args supplied")
        fmt.Println("arg error - usage: ./minimize <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[1]
    refGenome := parseGenome(genomeName)
    fmt.Println("number of chromosomes parsed:", len(refGenome.Chroms))

    /*
    fmt.Println()
    for k, v := range refGenome.Chroms {
        for k_, v_ := range v {
            fmt.Printf("refGenome.Chroms[%s][%s] = %s . . . %s\n", k, k_, v_[:10], v_[len(v_)-10:])
            fmt.Printf("len(refGenome.Chroms[%s][%s]) = %d\n", k, k_, len(v_))
        }
        fmt.Println()
    }
    */

    rawMatchesBytes, err := ioutil.ReadFile("dm3/dm3.fa.out")
    checkError(err)
    rawMatches := string(rawMatchesBytes)
    numLines, matchLines := lines(rawMatches)
    // drop header
    fmt.Println("number of parsed matchLines (including header):", numLines)
    matchLines = matchLines[3:]
    fmt.Println("number of parsed matchLines (excluding header):", len(matchLines))
    matches := parseMatches(matchLines)
    fmt.Println("number of parsed matches:", len(matches))

    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    var repeats []Repeat
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    var repeatsSize int64 = 1
    for i := range matches {
        repeatsSize = max(repeatsSize, matches[i].RepeatID + 1)
    }
    repeats = make([]Repeat, repeatsSize)
    // the real repeats begin at 1, we store root at 0
    repeats[0] = Repeat{0, []string{"root"}, "root"}
    //root := &repeats[0]

    // maps a repeat's category to its ID
    repeatMap := make(map[string](int64))
    // we now assign the actual repeats
    for i := 1; i < len(matches); i++ {
        id := matches[i].RepeatID
        // don't bother overwriting
        if repeats[id].ID == 0 {
            repeats[id].ID = id
            repeats[id].Class = append(strings.Split(matches[i].RepeatClass, "/"), matches[i].RepeatName)
            // "Other" or "Unknown" are meaningless categories, so all its children should be connected to root
            // !!! we should probably check in the future if "Other" or "Unknown" appear at ancestry levels other than the first
            // that all repeats descend from root is implicit - root is excluded from the slice
            if repeats[id].Class[0] == "Other" || repeats[id].Class[0] == "Unknown" {
                repeats[id].Class = repeats[id].Class[1:]
            } else {
                repeats[id].Class = repeats[id].Class
            }
            repeats[id].FullName = strings.Join([]string{matches[i].RepeatClass, matches[i].RepeatName}, "/")
        }
        repeatMap[strings.Join(repeats[id].Class, "/")] = id
    }

    classRoot, classNodes := GetClassTree(repeats)

    classRoot.PrintTree()
    fmt.Println("number of ClassNodes:", len(classNodes))

    fmt.Println("GetLCA(classNodes, classNodes['DNA/P/Galileo_DM'], classNodes['DNA/TcMar-Pogo/POGO']):", GetLCA(classNodes, classNodes["DNA/P/Galileo_DM"], classNodes["DNA/TcMar-Pogo/POGO"]))

    // matchBares := parseMatchBares(matchLines)

    fmt.Println("getMinimizer(\"ataggatcacgac\", 4) =", getMinimizer("ataggatcacgac", 4))
    fmt.Println("reverseComplement(\"aaAtGctACggT\") =", reverseComplement("aaAtGctACggT"))
}
