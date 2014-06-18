/*
    WARNING!!! This program is currently under development and may be buggy or broken.

    Must change kmer map index from kmer to minimizer.

    Still needs full support for below data spec.

    Int sizes should be reviewed for memory efficiency.

    A barebones (at the moment) Go script for parsing and minimizing repeats

    The sole command line argument is the name of the reference genome (e.g. "dm3").

    This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
        * a RepeatMasker library containing:
            - the match library (e.g. "dm3.fa.out")
            - the alignment information (e.g. "dm3.fa.align")
        * one or more reference genome files in FASTA format with the suffix ".fa"
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
    RepeatClass []string
    RepeatStart int64
    RepeatEnd int64
    RepeatRemains int64
    RepeatID int64

    // not in parsed data file
    FullName string
}


/*
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
*/


type RefGenome struct {
    Name string
    // maps a chromosome name to a map of its sequences
    Chroms map[string](map[string]string)
}


type ClassTree struct {
    // maps all class names to a pointer to their struct
    // chosen because of this foible in golang: https://code.google.com/p/go/issues/detail?id=3117
    // may be changed in the future
    ClassNodes map[string](*ClassNode)
    // a pointer to the the class tree's root, used for recursive descent etc.
    Root *ClassNode
}


type Kmer struct {
    // should likely be renamed "Seq" or "KmerSeq"
    Kmer string
    MinOffset uint8
    IsRevComp bool
    LCA *ClassNode
    Count map[int64]int32
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


// returns the number of lines and a slice of the lines
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


func ParseMatches(genomeName string) []Match {
    // !!! Below string literal is a temporary solution!
    rawMatchesBytes, err := ioutil.ReadFile("dm3/dm3.fa.out")
    checkError(err)
    rawMatches := string(rawMatchesBytes)
    _, matchLines := lines(rawMatches)
    // drop header
    matchLines = matchLines[3:]
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

        // everything in this block is just vanilla trimming, converting, and error checking
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
        match.RepeatClass = append(strings.Split(strings.TrimSpace(rawVals[10]), "/"), strings.TrimSpace(rawVals[9]))
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

        // decrement match.SeqStart and match.RepeatStart so that they work from a start index of 0 rather than 1
        // that way, we can use them without modification in slices
        match.SeqStart--
        match.RepeatStart--


        // "Other" and "Unknown" classes are heirarchically meaningless and really just mean "root", so we remove them
        if match.RepeatClass[0] == "Other" || match.RepeatClass[0] == "Unknown" {
            match.RepeatClass = match.RepeatClass[1:]
        }
        match.FullName = strings.Join(match.RepeatClass, "/")

        matches = append(matches, match)
    }

    return matches
}


/*
// !!! WARNING - this function is currently obsolete
// it may be used later for spacial locality
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
*/


func ParseGenome(genomeName string) RefGenome {
    // most fields will be assigned separately
    refGenome := RefGenome{}
    refGenome.Name = genomeName

    chromFileInfos, err := ioutil.ReadDir(genomeName)
    checkError(err)
    // used below to store the two keys for RefGenome.Chroms
    refGenome.Chroms = make(map[string](map[string]string))
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


func GetRepeats(matches []Match) []Repeat {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    var repeatsSize int64 = 1
    for i := range matches {
        repeatsSize = max(repeatsSize, matches[i].RepeatID + 1)
    }
    repeats := make([]Repeat, repeatsSize)

    // maps a repeat's category to its ID
    repeatMap := make(map[string](int64))
    // we now assign the actual repeats
    for i := 0; i < len(matches); i++ {
        id := matches[i].RepeatID
        // don't bother overwriting
        if repeats[id].ID == 0 {
            repeats[id].ID = id
            repeats[id].Class = matches[i].RepeatClass
            // "Other" or "Unknown" are meaningless categories, so all its children should be connected to root
            // !!! we should probably check in the future if "Other" or "Unknown" appear at ancestry levels other than the first
            // that all repeats descend from root is implicit - root is excluded from the slice
            //if repeats[id].Class[0] == "Other" || repeats[id].Class[0] == "Unknown" {
            //    repeats[id].Class = repeats[id].Class[1:]
            //}
            repeats[id].FullName = matches[i].FullName
            repeatMap[repeats[id].FullName] = id
        }
    }
    return repeats
}


func getMinimizer(kmer string, m uint8) (uint8, bool) {
    var minOffset uint8 = 0
    isRevComp := false
    var currMin, possMin string
    var i uint8
    for i = 0; i <= uint8(len(kmer)) - m; i++ {
        currMin = kmer[minOffset:minOffset+m]
        possMin = kmer[i:i+m]
        if seqCompare(currMin, possMin) == 1 {
            minOffset = i
            isRevComp = false
        }
        if seqCompare(currMin, revComp(possMin)) == 1 {
            minOffset = i
            isRevComp = true
        }
    }
    return minOffset, isRevComp
}


func revComp(seq string) string {
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


func (classTree ClassTree) PrintTree() {
    classTree.Root.printTreeRec(0)
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
func GetClassTree(repeats []Repeat) ClassTree {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    classNodes := make(map[string](*ClassNode))
    // would be prettier if expanded
    root := &ClassNode{}
    classNodes["root"] = root
    root.Name = "root"
    // all but Name is left nil
    for i := 1; i < len(repeats); i++ {
        // ignore the null indices
        if repeats[i].ID != 0 {
            if repeats[i].FullName == "root" {
                fmt.Println(repeats[i])
            }
            // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
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
                    // first case handles primary classes, as root is implicit and not listed in thisClass
                    if j == 1 {
                        classNodes[thisClassName].Parent = root
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
    return ClassTree{classNodes, root}
}


func (classTree ClassTree) getLCA(cnA, cnB *ClassNode) *ClassNode {
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
    return classTree.ClassNodes["root"]
}


// The logic for determining the minimizer
// Currently, it uses simple lexicographic ordering
// It returns -1 if a < b, 0 if a == b, and 1 if a > b
// It also assumes that the two sequences are of equal length, ignoring any hanging bases. In the future, maybe an error should be reported for this.
func seqCompare(a, b string) int8 {
    for i := 0; i < min(len(a), len(b)); i++ {
        if a[i] < b[i] {
            return -1
        }
        if a[i] > b[i] {
            return 1
        }
    }
    return 0
}


func (refGenome RefGenome) PrintChromInfo() {
    fmt.Println()
    for k, v := range refGenome.Chroms {
        for k_, v_ := range v {
            fmt.Printf("refGenome.Chroms[%s][%s] = %s . . . %s\n", k, k_, v_[:10], v_[len(v_)-10:])
            fmt.Printf("len(refGenome.Chroms[%s][%s]) = %d\n", k, k_, len(v_))
        }
        fmt.Println()
    }
}


// returns a pointer to the supplied kmer with the supplied minimizer in the supplied map
// if it does not exist, nil is returned
func getKmer(kmers map[string]([]*Kmer), minimizer string, kmer string) *Kmer {
    for i := range kmers[minimizer] {
        if kmers[minimizer][i].Kmer == kmer {
            return kmers[minimizer][i]
        }
    }
    return nil
}


func Minimize(refGenome RefGenome, matches []Match, classTree ClassTree, k uint8, m uint8) map[string]([]*Kmer) {
    // maps to *Kmer, as we did with ClassNodes, because of this golang foible: https://code.google.com/p/go/issues/detail?id=3117
    kmers := make(map[string]([]*Kmer))
    // because it's a uint8, k can be at most 255 (the same, of course, goes for m)
    var start, end int64
    // minOffset holds the offset of the current minimizer
    var x, minOffset uint8
    var seq, currMin, possMin, kmer string
    var isRevComp bool
    var kmerStruct *Kmer
    for i := range matches {
        start, end = matches[i].SeqStart, matches[i].SeqEnd
        // !!! the reference genome is two-dimensional, but RepeatMasker only supplies one sequence name
        // we resolve this ambiguity with the assumption that each chromosome file contains only one sequence
        // this holds true, at least in dm3
        seq = refGenome.Chroms[matches[i].SeqName][matches[i].SeqName][start:end]
        for j := 0; j <= len(seq) - int(k); j++ {
            kmer = seq[j:j+int(k)]
            // fmt.Printf("kmer: seq[%d:%d] = %s\n", j, j+int(k), kmer)    // DELETE ME
            // we have to calculate the first minimizer from scratch
            minOffset, isRevComp = getMinimizer(kmer, m)
            // x is the start index of the current possMin, the minimizer we're testing on this loop
            for x = 0; x <= k - m; x++ {
                // below is a potential performance sink - it can be optimized out later if necessary
                // fmt.Printf("possMin: kmer[%d:%d]\n", x, (x+m))    // DELETE ME
                possMin = kmer[x:x+m]
                // fmt.Println("minOffset:", minOffset)    // DELETE ME
                // holds the current minimizer
                currMin = kmer[minOffset:minOffset+m]
                // in some cases we can use most of the calculations from the previous kmer
                // when the previous minimizer isn't in the current kmer, though, we have to start from scratch
                if minOffset < x {
                    minOffset, isRevComp = getMinimizer(kmer, m)
                } else {
                    // otherwise there are only two options to overtake the previous minimizer: the last possible minimizer, and its reverse complement
                    if seqCompare(currMin, possMin) == 1 {
                        minOffset = x
                        isRevComp = false
                    }
                    if seqCompare(currMin, revComp(possMin)) == 1 {
                        minOffset = x
                        isRevComp = true
                    }
                }
                // if the kmers already there, we just update its LCA and increment its counter
                currMin = kmer[minOffset:minOffset+m]
                kmerStruct = getKmer(kmers, currMin, kmer)
                if kmerStruct != nil {
                    if classTree.ClassNodes[matches[i].FullName] == nil {
                        fmt.Println("match has nil ClassNode:", matches[i].FullName)
                        for y := 0; y < len(matches[i].RepeatClass); y++ {
                            fmt.Printf("location of ClassNodes['%s']: %p\n", strings.Join(matches[i].RepeatClass[0:y+1], ""), classTree.ClassNodes[strings.Join(matches[i].RepeatClass[0:y+1], "")])
                        }
                    }
                    // this shows a bit of a wart in the data structures - do we want each match to store a pointer to its ClassNode?
                    kmerStruct.LCA = classTree.getLCA(kmerStruct.LCA, classTree.ClassNodes[matches[i].FullName])
                    kmerStruct.Count[matches[i].RepeatID]++
                } else {
                    kmers[currMin] = append(kmers[currMin], &Kmer{kmer,
                                                                  minOffset,
                                                                  isRevComp,
                                                                  classTree.ClassNodes[matches[i].FullName],
                                                                  map[int64]int32{matches[i].RepeatID: 1}})
                }
            }
        }
    }
    return kmers
}


func main() {

    if len(os.Args) != 2 {
        fmt.Println(len(os.Args), "args supplied")
        fmt.Println("arg error - usage: ./minimize <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[1]
    refGenome := ParseGenome(genomeName)
    matches := ParseMatches(genomeName)
    repeats := GetRepeats(matches)
    classTree := GetClassTree(repeats)
    minimizers := Minimize(refGenome, matches, classTree, 31, 15)
    fmt.Println("number of kmers minimized:", len(minimizers))

    // below are testing statements

    fmt.Println("number of chromosomes parsed:", len(refGenome.Chroms))

    classTree.PrintTree()
    fmt.Println("number of ClassNodes:", len(classTree.ClassNodes))

    //fmt.Println("getLCA(classTree.ClassNodes, classTree.ClassNodes['DNA/P/Galileo_DM'], classTree.ClassNodes['DNA/TcMar-Pogo/POGO']):", getLCA(classTree.ClassNodes, classTree.ClassNodes["DNA/P/Galileo_DM"], classTree.ClassNodes["DNA/TcMar-Pogo/POGO"]))

    fmt.Println("min(5, 7):", min(5, 7))
    fmt.Println("max(int64(5), int64(7)):", max(int64(5), int64(7)))
    // matchBares := parseMatchBares(matchLines)

    //fmt.Println("getMinimizer(\"ataggatcacgac\", 4) =", getMinimizer("ataggatcacgac", 4))
    fmt.Println("revComp(\"aaAtGctACggT\") =", revComp("aaAtGctACggT"))
}
