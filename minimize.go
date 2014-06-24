/*
   WARNING!!! This program is currently under development and may be buggy or broken.

   A barebones (at the moment) Go script for parsing and minimizing repeats

   This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
       * a RepeatMasker library containing:
           - the match library (e.g. "dm3.fa.out")
           - the alignment information (e.g. "dm3.fa.align")
       * one or more reference genome files in FASTA format with the suffix ".fa"

   Premature commenting is the root of all evil, and I have sinned. Please read comments skeptically.

   There is no need to be referencing and dereferencing maps across funcs and structs because they're reference types

   The IsRevComp values in the output is suspect; should probably double-check that that's working.

   Should reconsider what is a pointer and what is directly referenced

   Should probably add a warning for actual 2D ref genomes

   Should review where []byte is used and where string is used

   For caching efficiency, we should change the minimizer data structure to a map-indexed 1D slice of Kmers (not *Kmers). (This technique originated in Kraken.)

   Int sizes should be reviewed for memory efficiency.

   The sole command line argument is the name of the reference genome (e.g. "dm3").
*/

package main

import (
    "bufio"
    "bytes"
    "flag"
    "fmt"
    "io/ioutil"
    "log"
    "os"
    "runtime"
    "runtime/pprof"
    "sort"
    "strconv"
    "strings"
    //"time"
    //"math"
)

var err error

// uses 64-bit values, which is probably unnecessary for some fields
// this is done for simplicity (and strconv compatibility)
// it should probably be changed eventually
type Match struct {
    SW_Score      int64
    PercDiv       float64
    PercDel       float64
    PercIns       float64
    SeqName       string
    SeqStart      int64
    SeqEnd        int64
    SeqRemains    int64
    IsRevComp     bool
    RepeatClass   []string
    RepeatStart   int64
    RepeatEnd     int64
    RepeatRemains int64
    RepeatID      int64

    // not in parsed data file
    FullName string
}

type RepeatGenome struct {
    Name string
    // maps a chromosome name to a map of its sequences
    Chroms    map[string](map[string][]byte)
    Kmers []Kmer
    OffsetsToMin map[string]int
    Matches   Matches
    ClassTree ClassTree
    Repeats   Repeats
    RepeatMap map[string]int64
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
    Kmer      []byte
    Min       string
    MinOffset uint8
    M         uint8
    IsRevComp bool
    LCA       *ClassNode
    Count     map[int64]int32
}

type Repeat struct {
    ID       int64
    Class    []string
    FullName string
}

type ClassNode struct {
    Name     string
    Class    []string
    Parent   *ClassNode
    Children []*ClassNode
    IsRoot   bool
}

// type synonym, necessary to implement interfaces (e.g. sort)
type Kmers []*Kmer
type Minimizers map[string]Kmers

type Repeats []Repeat
type Matches []Match

//type Chroms map[string](map[string][]byte)

// probably not the best name option
type MinPair struct {
    Kmer  *Kmer
    Match *Match
}

func checkError(err error) {
    if err != nil {
        log.Fatal(err)
    }
}

// returns the number of lines and a slice of the lines
func lines(byteSlice []byte) (numLines int, lines [][]byte) {
    numLines = 0
    for i := range byteSlice {
        if byteSlice[i] == '\n' {
            numLines++
        }
    }
    lines = bytes.Split(byteSlice, []byte{'\n'})
    // drop the trailing newline's line if it's there
    if len(bytes.TrimSpace(lines[len(lines)-1])) == 0 {
        lines = lines[:len(lines)-1]
    }
    return numLines, lines
}

func parseMatches(genomeName string) Matches {
    filepath := strings.Join([]string{genomeName, "/", genomeName, ".fa.out"}, "")
    rawMatchesBytes, err := ioutil.ReadFile(filepath)
    checkError(err)
    _, matchLines := lines(rawMatchesBytes)
    // drop header
    matchLines = matchLines[3:]
    var matches Matches

    for i := range matchLines {
        matchLine := string(matchLines[i])
        rawVals := strings.Fields(matchLine)
        var match Match
        match.IsRevComp = rawVals[8] == "C"

        // remove enclosing parentheses
        // !!! in the future, checkes to ensure that the parentheses exist should be added
        // !!! it would also be sensible to check that rawVals[8] is either "C" or "+"
        rawVals[7] = rawVals[7][1 : len(rawVals[7])-1]
        if match.IsRevComp {
            rawVals[11] = rawVals[11][1 : len(rawVals[11])-1]
        } else {
            rawVals[13] = rawVals[13][1 : len(rawVals[13])-1]
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
                match.RepeatRemains, match.RepeatStart, match.RepeatRemains+(match.RepeatEnd-match.RepeatStart)
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

func parseGenome(genomeName string) map[string](map[string][]byte) {
    chromFileInfos, err := ioutil.ReadDir(genomeName)
    checkError(err)
    chroms := make(map[string](map[string][]byte))
    // used below to store the two keys for RepeatGenome.chroms
    for i := range chromFileInfos {
        chromFilename := strings.Join([]string{genomeName, chromFileInfos[i].Name()}, "/")
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align)
        if strings.HasSuffix(chromFilename, ".fa") {
            rawSeqBytes, err := ioutil.ReadFile(chromFilename)
            checkError(err)
            numLines, seqLines := lines(rawSeqBytes)

            // maps each sequence name in this chrom to a slice of its sequence's lines
            // the list is concatenated at the end for efficiency's sake
            seqMap := make(map[string][][]byte)
            var seqName string
            for i := 0; i < numLines; i++ {
                seqLine := bytes.TrimSpace(seqLines[i])
                if seqLine[0] == byte('>') {
                    seqName = string(bytes.TrimSpace(seqLine[1:]))
                } else {
                    seqMap[seqName] = append(seqMap[seqName], seqLine)
                }
            }
            // finally, we insert this map into the full map
            chromName := chromFilename[len(genomeName)+1 : len(chromFilename)-3]
            chroms[chromName] = make(map[string][]byte)
            for k, v := range seqMap {
                chroms[chromName][k] = bytes.ToLower(bytes.Join(v, []byte{}))
            }
        }
    }
    return chroms
}

func GenerateRepeatGenome(genomeName string, k, m uint8) *RepeatGenome {
    repeatGenome := new(RepeatGenome)
    repeatGenome.Name = genomeName
    repeatGenome.Chroms = parseGenome(genomeName)
    repeatGenome.Matches = parseMatches(genomeName)
    repeatGenome.Repeats, repeatGenome.RepeatMap = repeatGenome.Matches.getRepeats()
    repeatGenome.ClassTree = repeatGenome.Repeats.getClassTree()

    runtime.GOMAXPROCS(runtime.NumCPU())
    minimizers := MinimizeServer(repeatGenome, k, m, runtime.NumCPU())
    minimizers.Write(genomeName)

    for _, kmers := range *minimizers {
        sort.Sort(kmers)
    }
    repeatGenome.OffsetsToMin = make(map[string]int)
    minSeqs := []string{}
    for key, _ := range *minimizers {
        minSeqs = append(minSeqs, key)
    }
    for i := range minSeqs {
        repeatGenome.OffsetsToMin[minSeqs[i]] = len(repeatGenome.Kmers)
        for j := range (*minimizers)[minSeqs[i]] {
            repeatGenome.Kmers = append(repeatGenome.Kmers, *(*minimizers)[minSeqs[i]][j])
        }
    }
    return repeatGenome
}

func (matches Matches) getRepeats() (Repeats, map[string]int64) {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    var repeatsSize int64 = 1
    for i := range matches {
        repeatsSize = max64(repeatsSize, matches[i].RepeatID+1)
    }
    repeats := make(Repeats, repeatsSize)

    // maps a repeat's category to its ID
    repeatMap := make(map[string](int64))
    // we now assign the actual repeats
    for i := range matches {
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
    return repeats, repeatMap
}

func getMinimizer(kmer []byte, m uint8) (uint8, bool) {
    var possMin []byte
    var i uint8
    var k uint8 = uint8(len(kmer))
    currMin := kmer[:m]
    var minOffset uint8 = 0
    isRevComp := false
    for i = 0; i <= k-m; i++ {
        // below is a potential performance sink - it can be optimized out later if necessary
        possMin = kmer[i : i+m]
        // holds the current minimizer
        // in some cases we can use most of the calculations from the previous kmer
        // when the previous minimizer isn't in the current kmer, though, we have to start from scratch
        if isRevComp && !testRevComp(currMin, possMin) {
            minOffset = i
            isRevComp = false
            currMin = kmer[minOffset : minOffset+m]
        }
        if isRevComp && testTwoRevComps(possMin, currMin) {
            minOffset = i
            isRevComp = true
            currMin = kmer[minOffset : minOffset+m]
        }
        if !isRevComp && seqLessThan(possMin, currMin) {
            minOffset = i
            isRevComp = false
            currMin = kmer[minOffset : minOffset+m]
        }
        if !isRevComp && testRevComp(possMin, currMin) {
            minOffset = i
            isRevComp = true
            currMin = kmer[minOffset : minOffset+m]
        }
    }
    return minOffset, isRevComp
}

func revComp(seq []byte) []byte {
    var revCompSeq = make([]byte, 0, len(seq))
    for i := range seq {
        switch seq[len(seq)-i-1] {
        case 'a':
            revCompSeq = append(revCompSeq, 't')
            break
        case 't':
            revCompSeq = append(revCompSeq, 'a')
            break
        case 'c':
            revCompSeq = append(revCompSeq, 'g')
            break
        case 'g':
            revCompSeq = append(revCompSeq, 'c')
            break
        case 'n':
            revCompSeq = append(revCompSeq, 'n')
            break
        }
    }
    return revCompSeq
}

func max(a int, b int) int {
    if a > b {
        return a
    } else {
        return b
    }
}

func max64(a int64, b int64) int64 {
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
    classTree.Root.printTreeRec(0, true)
}

// doesn't print leaves
// prevents the terminal from being flooded with Unknowns, Others, and Simple Repeats
func (classTree ClassTree) PrintBranches() {
    classTree.Root.printTreeRec(0, false)
}

func (classNode *ClassNode) printTreeRec(indent int, printLeaves bool) {
    for i := 0; i < indent; i++ {
        fmt.Printf("\t")
    }
    fmt.Println(classNode.Class[len(classNode.Class)-1])
    for i := range classNode.Children {
        if printLeaves || len(classNode.Children[i].Children) > 0 {
            classNode.Children[i].printTreeRec(indent+1, printLeaves)
        }
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
func (repeats Repeats) getClassTree() ClassTree {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    classNodes := make(map[string](*ClassNode))
    // would be prettier if expanded
    root := &ClassNode{"root", []string{"root"}, nil, nil, true}
    classNodes["root"] = root
    // all but Name is left nil
    for i := 1; i < len(repeats); i++ {
        // ignore the null indices
        if repeats[i].ID != 0 {
            // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
            for j := 1; j <= len(repeats[i].Class); j++ {
                thisClass := repeats[i].Class[:j]
                thisClassName := strings.Join(thisClass, "/")
                var keyExists bool
                _, keyExists = classNodes[thisClassName]
                if !keyExists {
                    thisClassNode := new(ClassNode)
                    thisClassNode.Name = thisClassName
                    thisClassNode.Class = thisClass
                    thisClassNode.IsRoot = false
                    classNodes[thisClassName] = thisClassNode
                    // first case handles primary classes, as root is implicit and not listed in thisClass
                    if j == 1 {
                        classNodes[thisClassName].Parent = root
                    } else {
                        classNodes[thisClassName].Parent = classNodes[strings.Join(thisClass[:len(thisClass)-1], "/")]
                    }
                    parent := classNodes[thisClassName].Parent
                    if parent.Children == nil {
                        parent.Children = make([]*ClassNode, 0)
                    }
                    parent.Children = append(parent.Children, classNodes[thisClassName])
                }
            }
        }
    }
    return ClassTree{classNodes, root}
}

func (cn *ClassNode) getAncestry() []*ClassNode {
    // condition to prevent nil dereference of Root.Parent
    if cn.IsRoot {
        return []*ClassNode{}
    }
    ancestry := []*ClassNode{cn.Parent}
    for !ancestry[len(ancestry)-1].IsRoot {
        ancestry = append(ancestry, ancestry[len(ancestry)-1].Parent)
    }
    return ancestry
}

func (classTree ClassTree) getLCA(cnA, cnB *ClassNode) *ClassNode {
    ancestryA, ancestryB := cnA.getAncestry(), cnB.getAncestry()
    for i := range ancestryA {
        if classSliceContains(ancestryB, ancestryA[i]) {
            return ancestryA[i]
        }
    }
    // necessary for compilation - Root should be in the ancestry paths
    return classTree.Root
}

// The logic for determining the minimizer
// Currently, it uses simple lexicographic ordering
func seqLessThan(a, b []byte) bool {
    size := min(len(a), len(b))
    for i := 0; i < size; i++ {
        if a[i] < b[i] {
            return true
        }
        if a[i] > b[i] {
            return false
        }
    }
    return false
}

func (refGenome RepeatGenome) PrintChromInfo() {
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
func (minimizers Minimizers) getKmer(minimizer string, kmer []byte) *Kmer {
    for i := range minimizers[minimizer] {
        if bytes.Equal(minimizers[minimizer][i].Kmer, kmer) {
            return minimizers[minimizer][i]
        }
    }
    return nil
}

// tests the reverse complement of seq against currMin (assumedly the current minimizer)
// true is returned only if the rev. comp. is lexicographically less than currMin
// this is used so that a new string doesn't have to be dynamically allocated for the rev. comp.
func testRevComp(a, b []byte) bool {
    var revCompChar byte
    for i := range a {
        switch a[len(a)-i-1] {
        case 'a':
            revCompChar = 't'
            break
        case 't':
            revCompChar = 'a'
            break
        case 'c':
            revCompChar = 'g'
            break
        case 'g':
            revCompChar = 'c'
            break
        case 'n':
            revCompChar = 'n'
            break
        }
        if revCompChar < b[i] {
            return true
        } else if revCompChar > b[i] {
            return false
        }
    }
    return false
}

// returns true only if the reverse complement of the first sequence is smaller than that of the second
func testTwoRevComps(a, b []byte) bool {
    var charA, charB byte
    for i := len(a) - 1; i >= 0; i-- {
        switch a[i] {
        case 'a':
            charA = 't'
            break
        case 't':
            charA = 'a'
            break
        case 'c':
            charA = 'g'
            break
        case 'g':
            charA = 'c'
            break
        case 'n':
            charA = 'n'
            break
        }

        switch b[i] {
        case 'a':
            charB = 't'
            break
        case 't':
            charB = 'a'
            break
        case 'c':
            charB = 'g'
            break
        case 'g':
            charB = 'c'
            break
        case 'n':
            charB = 'n'
            break
        }

        if charA < charB {
            return true
        } else if charA > charB {
            return false
        }
    }
    return false
}

// needed for sort.Interface
func (kmers Kmers) Len() int {
    return len(kmers)
}

func (kmers Kmers) Swap(i, j int) {
    kmers[i], kmers[j] = kmers[j], kmers[i]
}

func (kmers Kmers) Less(i, j int) bool {
    return seqLessThan(kmers[i].Kmer, kmers[j].Kmer)
}

func boolToInt(a bool) int {
    if a {
        return 1
    } else {
        return 0
    }
}

// !!! using []bytes instead of strings would probably make this faster
// additionally, using strconv and the writer's writing methods would probably be much faster
// profiling will determine whether these optimizations are worthwhile
func (minimizers Minimizers) Write(genomeName string) {
    filename := strings.Join([]string{genomeName, ".mins"}, "")
    outfile, err := os.Create(filename)
    checkError(err)
    defer outfile.Close()
    writer := bufio.NewWriter(outfile)
    defer writer.Flush()

    for thisMin, kmers := range minimizers {
        _, err = fmt.Fprintf(writer, ">%s\n", thisMin)
        checkError(err)
        for i := range kmers {
            _, err = fmt.Fprintf(writer, "\t%s %d %d\n", kmers[i].Kmer, kmers[i].MinOffset, boolToInt(kmers[i].IsRevComp))
            checkError(err)
            for repeatID, count := range kmers[i].Count {
                _, err = fmt.Fprintf(writer, "\t\t%d %d\n", repeatID, count)
                checkError(err)
            }
        }
    }
}

// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func MinimizeThread(repeatGenome *RepeatGenome, k, m uint8, matchStart, matchEnd int64, c chan MinPair) {
    //startTime := time.Now()
    var start, end int64
    var seqName string
    var seq, possMin, currMin, kmer, realMin []byte
    var minOffset, x uint8
    var isRevComp bool

    for i := matchStart; i < matchEnd; i++ {
        start, end, seqName = repeatGenome.Matches[i].SeqStart, repeatGenome.Matches[i].SeqEnd, repeatGenome.Matches[i].SeqName
        // !!! the reference genome is two-dimensional, but RepeatMasker only supplies one sequence name
        // we resolve this ambiguity with the assumption that each chromosome file contains only one sequence
        // this holds true, at least in dm3
        seq = repeatGenome.Chroms[seqName][seqName][start:end]
        for j := 0; j <= len(seq)-int(k); j++ {
            kmer = seq[j : j+int(k)]
            // we have to calculate the first minimizer from scratch
            minOffset, isRevComp = getMinimizer(kmer, m)
            currMin = kmer[minOffset : minOffset+m]
            // x is the start index of the current possMin, the minimizer we're testing on this loop
            for x = 0; x <= k-m; x++ {
                // below is a potential performance sink - it can be optimized out later if necessary
                possMin = kmer[x : x+m]
                // holds the current minimizer
                // in some cases we can use most of the calculations from the previous kmer
                // when the previous minimizer isn't in the current kmer, though, we have to start from scratch
                if minOffset < x {
                    minOffset, isRevComp = getMinimizer(kmer, m)
                    currMin = kmer[minOffset : minOffset+m]
                } else { 
                    if isRevComp && !testRevComp(currMin, possMin) {
                        minOffset = x
                        isRevComp = false
                        currMin = kmer[minOffset : minOffset+m]
                    }
                    if isRevComp && testTwoRevComps(possMin, currMin) {
                        minOffset = x
                        isRevComp = true
                        currMin = kmer[minOffset : minOffset+m]
                    }
                    if !isRevComp && seqLessThan(possMin, currMin) {
                        minOffset = x
                        isRevComp = false
                        currMin = kmer[minOffset : minOffset+m]
                    }
                    if !isRevComp && testRevComp(possMin, currMin) {
                        minOffset = x
                        isRevComp = true
                        currMin = kmer[minOffset : minOffset+m]
                    }
                }
            }
            // computationally expensive, so we only do it once at the end
            if isRevComp {
                realMin = revComp(currMin)
            } else {
                realMin = currMin
            }

            c <- MinPair{&Kmer{kmer,
                string(realMin),
                minOffset,
                m,
                isRevComp,
                repeatGenome.ClassTree.ClassNodes[repeatGenome.Matches[i].FullName],
                map[int64]int32{repeatGenome.Matches[i].RepeatID: 1}},
                &repeatGenome.Matches[i]}
        }
    }
    //fmt.Println("Thread life:", time.Since(startTime))
}

func MinimizeServer(repeatGenome *RepeatGenome, k, m uint8, numCPU int) *Minimizers {
    c := make(chan MinPair, 1000000)
    var i, mStart, mEnd int64
    for i := 0; i < numCPU; i++ {
        mStart = int64(i * len(repeatGenome.Matches) / numCPU)
        mEnd = int64((i + 1) * len(repeatGenome.Matches) / numCPU)
        go MinimizeThread(repeatGenome, k, m, mStart, mEnd, c)
    }

    var kmer, existingKmer *Kmer
    var match *Match
    var minPair MinPair
    minimizers := make(Minimizers)
    numKmers := repeatGenome.numKmers(k)
    fmt.Println(numKmers, "\t\tkmers to process")
    for i = 0; i < numKmers; i++ {
        if i%500000 == 0 {
            fmt.Println(i, "kmers processed")
        }
        minPair = <-c
        kmer, match = minPair.Kmer, minPair.Match
        existingKmer = minimizers.getKmer(kmer.Min, kmer.Kmer)
        if existingKmer != nil {
            // this shows a bit of a wart in the data structures - do we want each match to store a pointer to its ClassNode?
            existingKmer.Count[match.RepeatID]++
            existingKmer.LCA = repeatGenome.ClassTree.getLCA(existingKmer.LCA, kmer.LCA)
        } else {
            minimizers[kmer.Min] = append(minimizers[kmer.Min], kmer)
        }
    }
    return &minimizers
}

func (repeatGenome RepeatGenome) numKmers(k uint8) int64 {
    var k_ = int64(k)
    var matchSize int64
    var numMins int64 = 0
    for i := range repeatGenome.Matches {
        matchSize = repeatGenome.Matches[i].SeqEnd - repeatGenome.Matches[i].SeqStart
        if matchSize < k_ {

        } else {
            numMins += matchSize - k_ + 1
        }
    }
    return numMins
}

func (kmer Kmer) getMin() []byte {
    if kmer.IsRevComp {
        return revComp(kmer.Kmer[kmer.MinOffset : kmer.MinOffset+kmer.M])
    } else {
        return kmer.Kmer[kmer.MinOffset : kmer.MinOffset+kmer.M]
    }
}

func main() {

    var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

    flag.Parse()
    if *cpuprofile != "" {
        f, err := os.Create(*cpuprofile)
        if err != nil {
            log.Fatal(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    } else {
        fmt.Println("profiler off")
    }

    if len(os.Args) < 2 {
        fmt.Println(len(os.Args), "args supplied")
        fmt.Println("arg error - usage: ./minimize <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[2]
    var k uint8 = 31
    var m uint8 = 15
    repeatGenome := GenerateRepeatGenome(genomeName, k, m)

    // below are testing statements

    fmt.Println()
    for k, v := range repeatGenome.Chroms {
        for k_, v_ := range v {
            fmt.Printf("chrom: %s\tseq: %s\t%s...%s\n", k, k_, v_[:20], v_[len(v_)-20:])
        }
    }
    fmt.Println()

    fmt.Println()
    fmt.Println("number of chromosomes parsed:", len(repeatGenome.Chroms))
    fmt.Println()

    repeatGenome.ClassTree.PrintBranches()
    fmt.Println()
    fmt.Println("number of ClassNodes:", len(repeatGenome.ClassTree.ClassNodes))
    fmt.Println()
    fmt.Println("seqLessThan('aacct', 'aacct'):", seqLessThan([]byte("aacct"), []byte("aacct")))
    fmt.Println("seqLessThan('aacca', 'aaccg'):", seqLessThan([]byte("aacca"), []byte("aaccg")))
    fmt.Println("seqLessThan('aacct', 'aacca'):", seqLessThan([]byte("aacct"), []byte("aacca")))
    fmt.Println()
    fmt.Println("testRevComp('aacct', 'aacct'):", testRevComp([]byte("aacct"), []byte("aacct")))
    fmt.Println("testRevComp('aacca', 'aaccg'):", testRevComp([]byte("aacca"), []byte("aaccg")))
    fmt.Println("testRevComp('aacct', 'aacca'):", testRevComp([]byte("aacct"), []byte("aacca")))
    fmt.Println()
    fmt.Println("testTwoRevComps('aacct', 'aacct'):", testTwoRevComps([]byte("aacct"), []byte("aacct")))
    fmt.Println("testTwoRevComps('aacca', 'aaccg'):", testTwoRevComps([]byte("aacca"), []byte("aaccg")))
    fmt.Println("testTwoRevComps('aacct', 'aacca'):", testTwoRevComps([]byte("aacct"), []byte("aacca")))

    //fmt.Println("classTree.getLCA(classTree.ClassNodes['DNA/TcMar-Mariner'], classTree.ClassNodes['DNA/TcMar-Tc1']):", classTree.getLCA(classTree.ClassNodes["DNA/TcMar-Mariner"], classTree.ClassNodes["DNA/TcMar-Tc1"]))

    fmt.Println()
    fmt.Println("min(5, 7):", min(5, 7))
    fmt.Println("max64(int64(5), int64(7)):", max64(int64(5), int64(7)))
    fmt.Println()

    minOffset, isRevComp := getMinimizer([]byte("tgctcctgtcatgcatacgcaggtcatgcat"), 15)
    fmt.Println("getMinimizer('tgctcctgtcatgcatacgcaggtcatgcat', 15):", minOffset, isRevComp)
    //fmt.Println("testTwoRevComps('atttat', 'aattat'):", testTwoRevComps([]byte("atttat"), []byte("atttaa")))
    //fmt.Println("testRevComp('atttat', 'atttaa'):", testRevComp([]byte("atttat"), []byte("atttaa")))

    //fmt.Println("getMinimizer(\"ataggatcacgac\", 4) =", getMinimizer("ataggatcacgac", 4))
    fmt.Println("revComp(\"aaAtGctACggT\") =", revComp([]byte("aaAtGctACggT")))

    fmt.Println()
    fmt.Println("number of CPUs available:", runtime.NumCPU())

    /*
    fmt.Println()
    fmt.Println("expected number of kmers:", repeatGenome.numKmers(k))
    var numKmers int32 = 0
    for _, v := range minimizers {
        for i := range v {
            for _, count := range v[i].Count {
                numKmers += count
            }
        }
    }
    fmt.Println("number of kmers minimized:", numKmers)
    fmt.Println()
    */
}
