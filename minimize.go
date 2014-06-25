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

type Match struct {
    SW_Score      int32
    PercDiv       float64
    PercDel       float64
    PercIns       float64
    SeqName       string
    SeqStart      uint64
    SeqEnd        uint64
    SeqRemains    uint64
    IsRevComp     bool
    RepeatClass   []string
    // in weird cases, RepeatStart can be negative, so they must be signed
    RepeatStart   int64
    RepeatEnd     int64
    RepeatRemains int64
    RepeatID      int32

    // not in parsed data file
    FullName      string
}

type RepeatGenome struct {
    Name         string
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    Chroms       map[string](map[string][]byte)
    K            uint8
    M            uint8
    Kmers        []Kmer
    OffsetsToMin map[string]int
    Matches      Matches
    ClassTree    ClassTree
    Repeats      Repeats
    RepeatMap    map[string]int32
}

type ClassTree struct {
    // maps all class names to a pointer to their node struct
    // we must use pointers because of this foible in golang: https://code.google.com/p/go/issues/detail?id=3117
    // if we didn't use pointers, we'd have to completely reassign the struct when adding parents, children, etc.
    ClassNodes map[string](*ClassNode)
    // a pointer to the the class tree's root, used for recursive descent etc.
    // we explicitly create the root (because RepeatMatcher doesn't)
    Root       *ClassNode
}

type Kmer struct {
    Seq       []byte
    // storing the full minimizer here definitely isn't ideal
    // however, it pushes the significant work (namely memory allocation) of doing so to the MinimizerThreads, ultimately increasing efficiency
    // we could explore returning a struct containing the minimizer and Kmer, and discarding the prior after using it as an index
    Min       string
    MinOffset uint8
    M         uint8
    IsRevComp bool
    LCA       *ClassNode
    CountMap  map[int32]int32
}

type Repeat struct {
    // assigned by RepeatMasker, in simple incremented order starting from 1
    // they are therefore not compatible across genomes
    // we give root ID = 0
    ID        int32
    // a list containing the repeat's ancestry path, from top down
    // root is implicit, and is therefore excluded from the list
    ClassList []string
    FullName  string
}

type ClassNode struct {
    Name     string
    Class    []string
    Parent   *ClassNode
    Children []*ClassNode
    IsRoot   bool
}

// type synonyms, necessary to implement interfaces (e.g. sort) and methods
type Kmers      []*Kmer
type Minimizers map[string]Kmers
type Repeats    []Repeat
type Matches    []Match

//type Chroms map[string](map[string][]byte)

// return type of the parallel MinimizeThreads
// necessary to return both together over a channel
type KmerAndMatch struct {
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
    lastLine := lines[len(lines) - 1]
    if len(lastLine) == 1 && lastLine[0] == '\n' {
        lines = lines[:len(lines)-1]
    }
    return numLines, lines
}

func parseMatches(genomeName string) Matches {
    // "my_genome_name"  ->  "my_genome_name/my_genome_name.fa.out"
    filepath := strings.Join([]string{genomeName, "/", genomeName, ".fa.out"}, "")
    rawMatchesBytes, err := ioutil.ReadFile(filepath)
    checkError(err)
    _, matchLines := lines(rawMatchesBytes)
    // drop header

    matchLines = matchLines[3:]
    var matches Matches
    var sw_Score, repeatID int64

    for i := range matchLines {
        matchLine := string(matchLines[i])
        rawVals := strings.Fields(matchLine)
        if len(rawVals) < 15 {
            fmt.Println("FATAL ERROR: match line supplied to parseMatches() less than 15 fields long.")
        }
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
        sw_Score, err = strconv.ParseInt(rawVals[0], 10, 32)
        checkError(err)
        match.SW_Score = int32(sw_Score)
        match.PercDiv, err = strconv.ParseFloat(rawVals[1], 64)
        checkError(err)
        match.PercDel, err = strconv.ParseFloat(rawVals[2], 64)
        checkError(err)
        match.PercIns, err = strconv.ParseFloat(rawVals[3], 64)
        checkError(err)
        match.SeqName = strings.TrimSpace(rawVals[4])
        match.SeqStart, err = strconv.ParseUint(rawVals[5], 10, 64)
        checkError(err)
        match.SeqEnd, err = strconv.ParseUint(rawVals[6], 10, 64)
        checkError(err)
        match.SeqRemains, err = strconv.ParseUint(rawVals[7], 10, 64)
        checkError(err)
        // match.IsComplement, rawVals[8], moved above
        match.RepeatClass = append(strings.Split(strings.TrimSpace(rawVals[10]), "/"), strings.TrimSpace(rawVals[9]))
        match.RepeatStart, err = strconv.ParseInt(rawVals[11], 10, 64)
        checkError(err)
        match.RepeatEnd, err = strconv.ParseInt(rawVals[12], 10, 64)
        checkError(err)
        match.RepeatRemains, err = strconv.ParseInt(rawVals[13], 10, 64)
        checkError(err)
        repeatID, err = strconv.ParseInt(rawVals[14], 10, 32)
        checkError(err)
        match.RepeatID = int32(repeatID)

        // necessary swaps to convert reverse complement repeat indexes to positive-strand indexes
        if match.IsRevComp {
            match.RepeatStart = match.RepeatRemains
            match.RepeatEnd = match.RepeatStart
            match.RepeatRemains = match.RepeatRemains+(match.RepeatEnd-match.RepeatStart)
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
        // "my_genome_name", "my_chrom_name"  ->  "my_genome_name/my_chrom_name"
        chromFilename := strings.Join([]string{genomeName, chromFileInfos[i].Name()}, "/")
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align) or anything else
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
                    if string(seqName) != genomeName {
                        fmt.Println("WARNING: reference genome is two-dimensional, containing sequences not named after their chromosome.")
                        fmt.Println("Because RepeatMasker supplied only one-dimensional indexing, this may cause unexpected behavior or program failure.")
                    }
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
    // we popoulate the RepeatGenome mostly with helper functions
    // we should consider whether it makes more sense for them to alter the object directly, than to return their results
    repeatGenome := new(RepeatGenome)
    repeatGenome.Name = genomeName
    repeatGenome.Chroms = parseGenome(genomeName)
    repeatGenome.Matches = parseMatches(genomeName)
    repeatGenome.Repeats, repeatGenome.RepeatMap = repeatGenome.Matches.getRepeats()
    repeatGenome.ClassTree = repeatGenome.Repeats.getClassTree()
    repeatGenome.K = k
    repeatGenome.M = m

    // calling the parallel minimizer and writing the result
    repeatGenome.MinimizeServer(true)

    return repeatGenome
}

func (matches Matches) getRepeats() (Repeats, map[string]int32) {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    var repeatsSize int32 = 1
    for i := range matches {
        repeatsSize = max32(repeatsSize, matches[i].RepeatID+1)
    }
    repeats := make(Repeats, repeatsSize)

    // maps a repeat's category to its ID
    repeatMap := make(map[string](int32))
    // we now assign the actual repeats
    for i := range matches {
        id := matches[i].RepeatID
        // don't bother overwriting
        if repeats[id].ID == 0 {
            repeats[id].ID = id
            repeats[id].ClassList = matches[i].RepeatClass
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

func max32(a int32, b int32) int32 {
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
    for j := range repeat.ClassList {
        for j_ := 0; j_ < j; j_++ {
            fmt.Printf("\t")
        }
        fmt.Printf("%s\n", repeat.ClassList[j])
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
            for j := 1; j <= len(repeats[i].ClassList); j++ {
                thisClass := repeats[i].ClassList[:j]
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
        if bytes.Equal(minimizers[minimizer][i].Seq, kmer) {
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
    return seqLessThan(kmers[i].Seq, kmers[j].Seq)
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
            _, err = fmt.Fprintf(writer, "\t%s %d %d\n", kmers[i].Seq, kmers[i].MinOffset, boolToInt(kmers[i].IsRevComp))
            checkError(err)
            for repeatID, count := range kmers[i].CountMap {
                _, err = fmt.Fprintf(writer, "\t\t%d %d\n", repeatID, count)
                checkError(err)
            }
        }
    }
}

// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func MinimizeThread(repeatGenome *RepeatGenome, k, m uint8, matchStart, matchEnd uint64, c chan KmerAndMatch) {
    //startTime := time.Now()
    var start, end uint64
    var seqName string
    var seq, possMin, currMin, kmer, realMin []byte
    var minOffset, x uint8
    var y int
    var isRevComp bool
    n := byte('n')

    for i := matchStart; i < matchEnd; i++ {
        start, end, seqName = repeatGenome.Matches[i].SeqStart, repeatGenome.Matches[i].SeqEnd, repeatGenome.Matches[i].SeqName
        // !!! the reference genome is two-dimensional, but RepeatMasker only supplies one sequence name
        // we resolve this ambiguity with the assumption that each chromosome file contains only one sequence
        // this holds true, at least in dm3
        seq = repeatGenome.Chroms[seqName][seqName][start:end]
        for j := 0; j <= len(seq)-int(k); j++ {
            kmer = seq[j : j+int(k)]
            // ignore kmers containing n's
            for y = range kmer {
                if kmer[y] == n {
                    continue
                }
            }
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

            c <- KmerAndMatch{&Kmer{kmer,
                string(realMin),
                minOffset,
                m,
                isRevComp,
                repeatGenome.ClassTree.ClassNodes[repeatGenome.Matches[i].FullName],
                map[int32]int32{repeatGenome.Matches[i].RepeatID: 1}},
                &repeatGenome.Matches[i]}
        }
    }
    //fmt.Println("Thread life:", time.Since(startTime))
}

func (repeatGenome RepeatGenome) MinimizeServer(writeMins bool) {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    runtime.GOMAXPROCS(runtime.NumCPU())
    numCPU := runtime.NumCPU()
    c := make(chan KmerAndMatch, 1000000)
    var mStart, mEnd uint64
    for i := 0; i < numCPU; i++ {
        mStart = uint64(i * len(repeatGenome.Matches) / numCPU)
        mEnd = uint64((i + 1) * len(repeatGenome.Matches) / numCPU)
        go MinimizeThread(&repeatGenome, repeatGenome.K, repeatGenome.M, mStart, mEnd, c)
    }

    var kmer, existingKmer *Kmer
    var match *Match
    var kmerAndMatch KmerAndMatch
    minimizers := make(Minimizers)
    numKmers := repeatGenome.numKmers(repeatGenome.K)
    fmt.Println(numKmers, "\t\tkmers to process")
    var i uint64
    for i = 0; i < numKmers; i++ {
        if i%500000 == 0 {
            fmt.Println(i, "kmers processed")
        }
        kmerAndMatch = <-c
        kmer, match = kmerAndMatch.Kmer, kmerAndMatch.Match
        existingKmer = minimizers.getKmer(kmer.Min, kmer.Seq)
        if existingKmer != nil {
            // this shows a bit of a wart in the data structures - do we want each match to store a pointer to its ClassNode?
            existingKmer.CountMap[match.RepeatID]++
            existingKmer.LCA = repeatGenome.ClassTree.getLCA(existingKmer.LCA, kmer.LCA)
        } else {
            minimizers[kmer.Min] = append(minimizers[kmer.Min], kmer)
        }
    }

    if writeMins {
        minimizers.Write(repeatGenome.Name)
    }

    for _, kmers := range minimizers {
        sort.Sort(kmers)
    }

    repeatGenome.OffsetsToMin = make(map[string]int)
    minSeqs := []string{}
    for key, _ := range minimizers {
        minSeqs = append(minSeqs, key)
    }
    for i := range minSeqs {
        repeatGenome.OffsetsToMin[minSeqs[i]] = len(repeatGenome.Kmers)
        for j := range minimizers[minSeqs[i]] {
            repeatGenome.Kmers = append(repeatGenome.Kmers, *minimizers[minSeqs[i]][j])
        }
    }
}

func (repeatGenome RepeatGenome) numKmers(k uint8) uint64 {
    var k_ = uint64(k)
    var matchSize uint64
    var numMins uint64 = 0
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
        return revComp(kmer.Seq[kmer.MinOffset : kmer.MinOffset+kmer.M])
    } else {
        return kmer.Seq[kmer.MinOffset : kmer.MinOffset+kmer.M]
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
