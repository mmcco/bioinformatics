/*
   WARNING!!! This program is currently under development and may be buggy or broken.

   A barebones (at the moment) Go script for parsing and minimizing repeats

   https://gist.github.com/anonymous/aaabfb1b9a6e8b10c687

   This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
       * a RepeatMasker library containing:
           - the match library (e.g. "dm3.fa.out")
           - the alignment information (e.g. "dm3.fa.align")
       * one or more reference genome files in FASTA format with the suffix ".fa"

   Premature commenting is the root of all evil, and I have sinned. Please read comments skeptically.

   I should probably change some variable names, like repeatGenome, to less verbose variants, and use more aliases.

   seqToInt and revCompToInt need minor cleanup and a potential name-change.

   Should reconsider what is a pointer and what is directly referenced

   All sequences containing Ns are currently ignored.

   We should consider taking end minimizers once the code base is more mature.

   We should also review how to deal with m <= len(match) < k.

   For caching efficiency, we should change the minimizer data structure to a map-indexed 1D slice of Kmers (not *Kmers). (This technique originated in Kraken.)

   Int sizes should be reviewed for memory efficiency.

   The sole command line argument is the name of the reference genome (e.g. "dm3").
*/

package main

import (
    "io"
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
    "sync"
    "time"
)

var err error

var DEBUG      *bool
var CPUPROFILE *bool
var MEMPROFILE *bool

// Match.SW_Score - Smith-Waterman score, describing the likeness to the repeat reference sequence
// Match.PercDiv
// Match.PercDel
// Match.PercIns
// Match.SeqName -  the reference genome file this match came from (typically the chromosome)
// Match.SeqStart -  the starting index (inclusive) in the reference genome
// Match.SeqEnd -  the ending index (exclusive) in the reference genome
// Match.SeqRemains
// Match.IsRevComp -  the match may be for the complement of the reference sequence
// Match.RepeatClass - the repeat's full ancestry, including its repeat class and repeat name (which are listed separately in the RepeatMasker output file)
// Match.RepeatStart-  the starting index in the repeat consensus sequence
// Match.RepeatEnd -  the ending sequence (exclusive) in the repeat consensus sequence
// Match.RepeatRemains
// Match.RepeatID - a numerical ID for the repeat type (starts at 1)
//     Match.Fullname - not in parsed data file - simply repeatClass concatenated - used for map indexing
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
    RepeatID      uint32

    FullName      string
    ClassNode     *ClassNode
}

type RepeatGenome struct {
    Name         string
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    Chroms       map[string](map[string]string)
    K            uint8
    M            uint8
    Kmers        []Kmer
    OffsetsToMin map[uint64]uint64
    SortedMins   Uint64Slice
    Matches      Matches
    ClassTree    ClassTree
    Repeats      Repeats
    RepeatMap    map[string]uint32
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
    SeqInt uint64
    LCA *ClassNode
}

type Repeat struct {
    // assigned by RepeatMasker, in simple incremented order starting from 1
    // they are therefore not compatible across genomes
    // we give root ID = 0
    ID        uint32
    // a list containing the repeat's ancestry path, from top down
    // root is implicit, and is therefore excluded from the list
    ClassList []string
    ClassNode *ClassNode
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
type MinMap     map[uint64]Kmers
type Repeats    []Repeat
type Matches    []Match

//type Chroms map[string](map[string][]byte)

type ThreadResponse struct {
    KmerInt      uint64
    Relative  *ClassNode
}

type MinCache struct {
    sync.Mutex
    Cache map[uint64]uint64
}

func checkError(err error) {
    if err != nil {
        log.Fatal(err)
    }
}

func seqToInt(seq string) uint64 {
    if len(seq) < 1 || len(seq) > 31 {
        panic("seqToInt() can only int-ize sequences where 0 < length < 32")
    }
    var seqInt uint64 = 0
    for i := 0; i < len(seq); i++ {
        seqInt = seqInt << 2
        switch seq[i] {
        case 'a':
            break
        case 'c':
            seqInt |= 1
            break
        case 'g':
            seqInt |= 2
            break
        case 't':
            seqInt |= 3
            break
        default:
            panic("byte other than 'a', 'c', 'g', or 't' passed to seqToInt()")
        }
    }
    return seqInt
}

func revCompToInt(seq string) uint64 {
    if len(seq) < 1 || len(seq) > 31 {
        panic("revCompToInt() can only int-ize sequences where 0 < length < 32")
    }
    var seqInt uint64 = 0
    for i := range seq {
        seqInt = seqInt << 2
        switch seq[len(seq) - (i + 1)] {
        case 'a':
            seqInt |= 3
            break
        case 'c':
            seqInt |= 2
            break
        case 'g':
            seqInt |= 1
            break
        case 't':
            break
        default:
            panic("byte other than 'a', 'c', 'g', or 't' passed to revCompToInt()")
        }
    }
    return seqInt
}

func intRevComp(posInt uint64, seqLen uint8) uint64 {
    var revComp uint64 = 0
    var i uint8
    for i = 0; i < seqLen; i++ {
        // should execute before every loop but the first
        if i != 0 {
            revComp <<= 2
            posInt >>= 2
        }
        switch posInt & 3 {
        case 0:
            revComp |= 3
            break
        case 1:
            revComp |= 2
            break
        case 2:
            revComp |= 1
            break
        case 3:
            break
        default:
            panic("bit manipulation logic error in intRevComp()")
        }
    }
    return revComp
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
    // drop the trailing newlines
    newline := []byte("\n")
    for lastLine := lines[len(lines) - 1]; len(lines) > 0 && (len(lastLine) == 0 || bytes.Equal(lastLine, newline)); lastLine = lines[len(lines) - 1] {
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
    repeatMap := make(map[uint32]string)

    var maxRepeatID uint32 = 1
    for i := range matchLines {
        matchLine := string(matchLines[i])
        rawVals := strings.Fields(matchLine)
        if len(rawVals) < 15 {
            fmt.Printf("FATAL ERROR: match line supplied to parseMatches() less than 15 fields long (has %d fields and length %d):\n", len(rawVals), len(matchLine))
            log.Fatal(matchLine)
        }
        repeatID, err := strconv.ParseUint(strings.TrimSpace(rawVals[14]), 10, 32)
        checkError(err)
        maxRepeatID = maxU32(maxRepeatID, uint32(repeatID))
    }
    var customRepeatID uint32 = maxRepeatID + 1

    for i := range matchLines {
        matchLine := string(matchLines[i])
        rawVals := strings.Fields(matchLine)
        if len(rawVals) < 15 {
            fmt.Printf("FATAL ERROR: match line supplied to parseMatches() less than 15 fields long (has %d fields and length %d):\n", len(rawVals), len(matchLine))
            log.Fatal(matchLine)
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
        match.RepeatID = uint32(repeatID)

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

        if repeatMap[match.RepeatID] != match.FullName {
            match.RepeatID = customRepeatID
            customRepeatID++
        }

        matches = append(matches, match)
    }
    return matches
}

func parseGenome(genomeName string) map[string](map[string]string) {
    chromFileInfos, err := ioutil.ReadDir(genomeName)
    checkError(err)
    warned := false
    chroms := make(map[string](map[string]string))
    // used below to store the two keys for RepeatGenome.chroms
    for i := range chromFileInfos {
        // "my_genome_name", "my_chrom_name"  ->  "my_genome_name/my_chrom_name"
        chromFilename := chromFileInfos[i].Name()
        chromFilepath := strings.Join([]string{genomeName, chromFilename}, "/")
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align) or anything else
        if strings.HasSuffix(chromFilepath, ".fa") {
            rawSeqBytes, err := ioutil.ReadFile(chromFilepath)
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
                    if !warned && seqName != chromFilename[:len(chromFilename) - 3] {
                        fmt.Println("WARNING: reference genome is two-dimensional, containing sequences not named after their chromosome.")
                        fmt.Println("Because RepeatMasker supplied only one-dimensional indexing, this may cause unexpected behavior or program failure.")
                        fmt.Println("seqName:", seqName, "\tlen(seqName):", len(seqName))
                        fmt.Println("chrom name:", chromFilename[:len(chromFilename) - 3], "\tlen(chrom name):", len(chromFilename) - 3)
                        warned = true
                    }
                } else {
                    seqMap[seqName] = append(seqMap[seqName], seqLine)
                }
            }
            // finally, we insert this map into the full map
            chromName := chromFilepath[len(genomeName)+1 : len(chromFilepath)-3]
            chroms[chromName] = make(map[string]string)
            for k, v := range seqMap {
                chroms[chromName][k] = strings.ToLower(string(bytes.Join(v, []byte{})))
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
    repeatGenome.getRepeats()
    repeatGenome.getClassTree()
    repeatGenome.K = k
    repeatGenome.M = m

    // calling the parallel minimizer and writing the result
    repeatGenome.GetKrakenSlice(true)

    return repeatGenome
}

func (repeatGenome *RepeatGenome) getRepeats() {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    var repeatsSize uint32 = 1
    for i := range repeatGenome.Matches {
        repeatsSize = maxU32(repeatsSize, repeatGenome.Matches[i].RepeatID+1)
    }
    repeatGenome.Repeats = make(Repeats, repeatsSize)

    // maps a repeat's category to its ID
    repeatGenome.RepeatMap = make(map[string](uint32))
    // we now assign the actual repeats
    for i := range repeatGenome.Matches {
        id := repeatGenome.Matches[i].RepeatID
        // don't bother overwriting
        if repeatGenome.Repeats[id].ID == 0 {
            repeatGenome.Repeats[id].ID = id
            repeatGenome.Repeats[id].ClassList = repeatGenome.Matches[i].RepeatClass
            repeatGenome.Repeats[id].FullName = repeatGenome.Matches[i].FullName
            repeatGenome.RepeatMap[repeatGenome.Repeats[id].FullName] = id
        }
    }
}

func getMinimizer(kmer string, m uint8) (uint8, uint64) {
    var k uint8 = uint8(len(kmer))
    numPossMins := k - (m - 1)
    if numPossMins < 1 || m < 1 {
        panic("getMinimizer(): m must be <= len(k) and > 0")
    }

    // stores the index of the leftmost base included in the minimizer
    var currOffset uint8 = 0
    currMin := seqToInt(kmer[:m])
    possMin := currMin
    // i indexes the base to be added in a given iteration
    for i := m; i < k; i++ {
        // overflow off the first base
        possMin <<= 64 - 2 * (m - 1)
        // return to proper alignment
        possMin >>= 64 - 2 * m
        // update new char
        switch kmer[i] {
        case 'a':
            break
        case 'c':
            possMin |= 1
            break
        case 'g':
            possMin |= 2
            break
        case 't':
            possMin |= 3
            break
        default:
            panic("getMinimizer(): kmer supplied containing byte other than 'a', 'c', 'g', or 't'")
        }
        if possMin < currMin {
            currMin = possMin
            currOffset = uint8(i) - (m - 1)
        }
    }

    // now we test the reverse complements

    possMin = revCompToInt(kmer[k - m : ])
    if possMin < currMin {
        currMin = possMin
        currOffset = k - m
    }
    // we now work right-to-left
    // i indexes the char to be added in a given iteration
    // must use i of type int because uint8 causes overflow
    for i := int(k) - (int(m) + 1); i >= 0; i-- {
        // overflow off the first base
        possMin <<= 64 - 2 * (m - 1)
        // return to proper alignment
        possMin >>= 64 - 2 * m
        // update new char
        switch kmer[i] {
        case 'a':
            possMin |= 3
            break
        case 'c':
            possMin |= 2
            break
        case 'g':
            possMin |= 1
            break
        case 't':
            break
        default:
            panic("getMinimizer(): kmer supplied containing byte other than 'a', 'c', 'g', or 't'")
        }
        if possMin < currMin {
            currMin = possMin
            currOffset = uint8(i)
        }
    }
    return currOffset, currMin
}

func revComp(seq string) string {
    var revCompSeq = make([]byte, 0, len(seq))
    // !!! DO NOT use range notation here - doing so is inefficient, converting the string to a []rune
    for i := 0; i < len(seq); i++ {
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
    return string(revCompSeq)
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

func maxU32(a uint32, b uint32) uint32 {
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

func (repeat *Repeat) Print() {
    for j := range repeat.ClassList {
        for j_ := 0; j_ < j; j_++ {
            fmt.Printf("\t")
        }
        fmt.Printf("%s\n", repeat.ClassList[j])
    }
    fmt.Println()
}

func (classTree *ClassTree) PrintTree() {
    classTree.Root.printTreeRec(0, true)
}

// doesn't print leaves
// prevents the terminal from being flooded with Unknowns, Others, and Simple Repeats
func (classTree *ClassTree) PrintBranches() {
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

func (repeatGenome *RepeatGenome) getClassTree() {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    repeatGenome.ClassTree.ClassNodes = make(map[string](*ClassNode))
    // would be prettier if expanded
    repeatGenome.ClassTree.Root = &ClassNode{"root", []string{"root"}, nil, nil, true}
    repeatGenome.ClassTree.ClassNodes["root"] = repeatGenome.ClassTree.Root
    for i := 1; i < len(repeatGenome.Repeats); i++ {
        // ignore the null indices
        if repeatGenome.Repeats[i].ID != 0 {
            // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
            for j := 1; j <= len(repeatGenome.Repeats[i].ClassList); j++ {
                thisClass := repeatGenome.Repeats[i].ClassList[:j]
                thisClassName := strings.Join(thisClass, "/")
                _, keyExists := repeatGenome.ClassTree.ClassNodes[thisClassName]
                if !keyExists {
                    thisClassNode := new(ClassNode)
                    thisClassNode.Name = thisClassName
                    thisClassNode.Class = thisClass
                    thisClassNode.IsRoot = false
                    repeatGenome.ClassTree.ClassNodes[thisClassName] = thisClassNode
                    // first case handles primary classes, as root is implicit and not listed in thisClass
                    if j == 1 {
                        repeatGenome.ClassTree.ClassNodes[thisClassName].Parent = repeatGenome.ClassTree.Root
                    } else {
                        repeatGenome.ClassTree.ClassNodes[thisClassName].Parent = repeatGenome.ClassTree.ClassNodes[strings.Join(thisClass[:len(thisClass)-1], "/")]
                    }
                    parent := repeatGenome.ClassTree.ClassNodes[thisClassName].Parent
                    if parent.Children == nil {
                        parent.Children = make([]*ClassNode, 0)
                    }
                    parent.Children = append(parent.Children, repeatGenome.ClassTree.ClassNodes[thisClassName])
                }
            }
        }
    }

    for i := range repeatGenome.Repeats {
        repeatGenome.Repeats[i].ClassNode = repeatGenome.ClassTree.ClassNodes[repeatGenome.Repeats[i].FullName]
    }

    for i := range repeatGenome.Matches {
        repeatGenome.Matches[i].ClassNode = repeatGenome.ClassTree.ClassNodes[repeatGenome.Matches[i].FullName]
        if repeatGenome.Matches[i].ClassNode == nil {
            fmt.Println(repeatGenome.Matches[i].FullName)
            log.Fatal("nil match ClassNode")
        }
    }
}

// returns ancestry list, from self to root
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

func (classTree *ClassTree) getLCA(cnA, cnB *ClassNode) *ClassNode {
    ancestryA := cnA.getAncestry()
    cnBWalker := cnB
    var i int
    for cnBWalker != classTree.Root {
        for i = 0; i < len(ancestryA); i++ {
            if cnBWalker == ancestryA[i] {
                return cnBWalker
            }
        }
        cnBWalker = cnBWalker.Parent
    }
    // necessary for compilation - Root should be in the ancestry paths
    return classTree.Root
}

// The logic for determining the minimizer
// Currently, it uses simple lexicographic ordering
func seqLessThan(a, b string) bool {
    // min function manually inlined for speed - dubious
    var size int
    if len(a) < len(b) {
        size = len(a)
    } else {
        size = len(b)
    }
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

func (repeats *Repeats) Write(filename string) {
    outfile, err := os.Create(filename)
    checkError(err)
    defer outfile.Close()

    for i := range *repeats {
        if int((*repeats)[i].ID) == i {
            fmt.Fprintf(outfile, "%d %s\n", (*repeats)[i].ID, (*repeats)[i].FullName)
        }
    }
}

func (refGenome *RepeatGenome) PrintChromInfo() {
    fmt.Println()
    for k, v := range refGenome.Chroms {
        for k_, v_ := range v {
            fmt.Printf("refGenome.Chroms[%s][%s] = %s . . . %s\n", k, k_, v_[:10], v_[len(v_)-10:])
            fmt.Printf("len(refGenome.Chroms[%s][%s]) = %d\n", k, k_, len(v_))
        }
        fmt.Println()
    }
}

// needed for sort.Interface
func (kmers Kmers) Len() int {
    return len(kmers)
}

func (kmers Kmers) Swap(i, j int) {
    kmers[i], kmers[j] = kmers[j], kmers[i]
}

func (kmers Kmers) Less(i, j int) bool {
    return kmers[i].SeqInt < kmers[j].SeqInt
}

// needed for sort.Interface
type Uint64Slice []uint64

func (uint64s Uint64Slice) Len() int {
    return len(uint64s)
}

func (uint64s Uint64Slice) Swap(i, j int) {
    uint64s[i], uint64s[j] = uint64s[j], uint64s[i]
}

func (uint64s Uint64Slice) Less(i, j int) bool {
    return uint64s[i] < uint64s[j]
}
func boolToInt(a bool) int {
    if a {
        return 1
    } else {
        return 0
    }
}

// a saner way of doing this would be to allocate a single k-long []byte and have a function populate it before printing
func (repeatGenome RepeatGenome) WriteMins() error {
    k := repeatGenome.K
    m := repeatGenome.M
    filename := strings.Join([]string{repeatGenome.Name, ".mins"}, "")
    outfile, err := os.Create(filename)
    if err != nil {
        return err
    }
    defer outfile.Close()
    writer := bufio.NewWriter(outfile)
    defer writer.Flush()

    var kmer *Kmer
    minIndex := 0
    thisMin := repeatGenome.SortedMins[minIndex]
    thisMinsBoundary := repeatGenome.OffsetsToMin[thisMin]
    _, err = fmt.Fprint(writer, ">")
    if err != nil {
        return err
    }
    err = writeSeqInt(writer, thisMin, m)
    if err != nil {
        return err
    }
    _, err = fmt.Fprint(writer, "\n")
    if err != nil {
        return err
    }
    for i := range repeatGenome.Kmers {
        if uint64(i) > thisMinsBoundary {
            minIndex++
            thisMin = repeatGenome.SortedMins[minIndex]
            thisMinsBoundary = repeatGenome.OffsetsToMin[thisMin]
            _, err = fmt.Fprint(writer, ">")
            if err != nil {
                return err
            }
            err = writeSeqInt(writer, thisMin, m)
            if err != nil {
                return err
            }
            _, err = fmt.Fprint(writer, "\n")
            if err != nil {
                return err
            }
        }
        kmer = &repeatGenome.Kmers[i]
        _, err = fmt.Fprint(writer, "\t")
        if err != nil {
            return err
        }
        err = writeSeqInt(writer, kmer.SeqInt, k)
        if err != nil {
            return err
        }
        _, err = fmt.Fprintf(writer, " %s\n", kmer.LCA.Name)
        if err != nil {
            return err
        }
    }
    return nil
}

func printSeqInt(seqInt uint64, seqLen uint8) {
    var i uint8
    for i = 0; i < seqLen; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (seqLen - i - 1))) & 3 {
        case 0:
            fmt.Print("a")
            break
        case 1:
            fmt.Print("c")
            break
        case 2:
            fmt.Print("g")
            break
        case 3:
            fmt.Print("t")
            break
        default:
            panic("error in printSeqInt() base selection")
        }
    }
}

func writeSeqInt(writer io.ByteWriter, seqInt uint64, seqLen uint8) error {
    var i uint8
    for i = 0; i < seqLen; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (seqLen - i - 1))) & 3 {
        case 0:
            err = writer.WriteByte('a')
            break
        case 1:
            err = writer.WriteByte('c')
            break
        case 2:
            err = writer.WriteByte('g')
            break
        case 3:
            err = writer.WriteByte('t')
            break
        default:
            panic("error in printSeqInt() base selection")
        }
        if err != nil {
            return err
        }
    }
    return nil
}

// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func (repeatGenome *RepeatGenome) minimizeThread(minCache *MinCache, matchStart, matchEnd uint64, c chan ThreadResponse) {
    var startTime time.Time
    if *DEBUG {
        startTime = time.Now()
    }

    var match *Match
    var seqLen uint64
    k := repeatGenome.K
    k_ := uint64(k)
    m := repeatGenome.M
    //m_ := uint64(m)
    var currOffset, x uint8
    var seq, kmerSeq string
    var kmerInt, possMin, currMin uint64
    var exists bool

    for i := matchStart; i < matchEnd; i++ {
        match = &repeatGenome.Matches[i]
        seq = repeatGenome.Chroms[match.SeqName][match.SeqName]
        seqLen = match.SeqEnd - match.SeqStart
        // for now, we will ignore matches too short to be traditionally minimized
        if seqLen < k_ {
            continue
        }

        KmerLoop:
        for j := match.SeqStart; j <= match.SeqEnd - k_; j++ {
            kmerSeq = seq[j : j + k_]
            // we begin by skipping any kmers containing n's
            for x = 0; x < k; x++ {
                if kmerSeq[x] == byte('n') {
                    j += uint64(x)
                    continue KmerLoop
                }
            }
            kmerInt = seqToInt(kmerSeq)
            minCache.Lock()
            _, exists = minCache.Cache[kmerInt]
            minCache.Unlock()
            if exists {
                c <- ThreadResponse{kmerInt, match.ClassNode}
            } else if j == match.SeqStart || uint64(currOffset) < j {
                currOffset, currMin = getMinimizer(kmerSeq, m)
                minCache.Lock()
                minCache.Cache[kmerInt] = currMin
                minCache.Unlock()

                c <- ThreadResponse{kmerInt, match.ClassNode}
            } else {
                possMin = seqToInt(kmerSeq[k - m : ])
                if possMin < currMin {
                    currMin = possMin
                    currOffset = k - m
                }
                possMin = revCompToInt(kmerSeq[k - m : ])
                if possMin < currMin {
                    currMin = possMin
                    currOffset = k - m
                }
                minCache.Lock()
                minCache.Cache[kmerInt] = currMin
                minCache.Unlock()

                c <- ThreadResponse{kmerInt, match.ClassNode}
            }
        }
    }

    if *DEBUG {
        fmt.Println("Thread life:", time.Since(startTime))
    }
}

func (repeatGenome *RepeatGenome) GetKrakenSlice(writeMins bool) {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    numCPU := runtime.NumCPU()
    if *DEBUG {
        fmt.Printf("using %d CPUs\n", numCPU)
    }
    runtime.GOMAXPROCS(numCPU)
    c := make(chan ThreadResponse, 100000)
    var mStart, mEnd uint64
    //minMap := make(MinMap)
    kmerMap := make(map[uint64]*Kmer)
    minCache := new(MinCache)
    minCache.Cache = make(map[uint64]uint64)

    numKmers := repeatGenome.numKmers()
    fmt.Println("expecting", numKmers, "kmers")

    for i := 0; i < numCPU; i++ {
        mStart = uint64(i * len(repeatGenome.Matches) / numCPU)
        mEnd = uint64((i + 1) * len(repeatGenome.Matches) / numCPU)
        go repeatGenome.minimizeThread(minCache, mStart, mEnd, c)
    }

    // below is the atomic section of minimizing, which is parallel
    // this seems to be the rate-limiting section, as threads use only ~6-7 CPU-equivalents
    // it should therefore be optimized before other sections
    var kmer *Kmer
    var response ThreadResponse
    var i, kmerInt uint64
    var relative *ClassNode
    var exists bool

    for i = 0; i < numKmers; i++ {
        if i%5000000 == 0 {
            fmt.Println(i / 1000000, "million kmers processed")
        }

        response = <- c
        kmerInt, relative = response.KmerInt, response.Relative

        if relative == nil {
            panic("nil relative returned by thread")
        }

        kmer, exists = kmerMap[kmerInt]
        if exists {
            kmer.LCA = repeatGenome.ClassTree.getLCA(kmer.LCA, relative)
        } else {
            kmerMap[kmerInt] = &Kmer{kmerInt, relative}
        }
    }

    fmt.Println("all minimizers generated")

    var kmers Kmers
    minMap := make(map[uint64]Kmers)
    for kmerInt, minimizer := range minCache.Cache {
        kmers, exists = minMap[minimizer]
        if exists {
            kmers = append(kmers, kmerMap[kmerInt])
        } else {
            minMap[minimizer] = Kmers{kmerMap[kmerInt]}
        }
    }

    repeatGenome.SortedMins = make(Uint64Slice, 0, len(minMap))
    for minimizer, _ := range minMap {
        repeatGenome.SortedMins = append(repeatGenome.SortedMins, minimizer)
    }
    sort.Sort(repeatGenome.SortedMins)

    var currOffset uint64 = 0
    repeatGenome.OffsetsToMin = make(map[uint64]uint64)
    repeatGenome.Kmers = []Kmer{}
    for i := range repeatGenome.SortedMins {
        repeatGenome.OffsetsToMin[repeatGenome.SortedMins[i]] = currOffset
        currOffset += uint64(len(minMap[repeatGenome.SortedMins[i]]))
        for j := range minMap[repeatGenome.SortedMins[i]] {
            repeatGenome.Kmers = append(repeatGenome.Kmers, *minMap[repeatGenome.SortedMins[i]][j])
        }
    }

    if writeMins {
        err = repeatGenome.WriteMins()
        checkError(err)
    }

    if *MEMPROFILE {
        f, err := os.Create(repeatGenome.Name + ".memprof")
        checkError(err)
        pprof.WriteHeapProfile(f)
        f.Close()
    }
}

func (repeatGenome *RepeatGenome) numKmers() uint64 {
    var k = int(repeatGenome.K)
    var numKmers uint64 = 0
    var seq string
    seqs := []string{}
    var match *Match

    splitOnN := func(c rune) bool {return c == 'n'}

    for i := range repeatGenome.Matches {
        match = &repeatGenome.Matches[i]
        seq = repeatGenome.Chroms[match.SeqName][match.SeqName][match.SeqStart : match.SeqEnd]
        seqs = strings.FieldsFunc(seq, splitOnN)
        for j := range seqs {
            if len(seqs[j]) >= k {
                numKmers += uint64(len(seqs[j]) - k + 1)
            }
        }
    }
    return numKmers
}

func main() {

    if len(os.Args) < 2 {
        fmt.Println("arg error - usage: ./minimize <flags> <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[len(os.Args) - 1]

    CPUPROFILE = flag.Bool("cpuprof", false, "write cpu profile to file")
    MEMPROFILE = flag.Bool("memprof", false, "write memory profile to this file")
    DEBUG = flag.Bool("debug", false, "run and print debugging tests")
    k_arg := flag.Uint("k", 31, "kmer length")
    m_arg := flag.Uint("m", 15, "minimizer length")
    flag.Parse()

    if *CPUPROFILE {
        fmt.Println("CPU profiler enabled")
        f, err := os.Create(genomeName + ".cpuprof")
        checkError(err)
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    } else {
        fmt.Println("CPU profiler disabled")
    }

    if *MEMPROFILE {
        fmt.Println("memory profiler enabled")
    } else {
        fmt.Println("memory profiler disabled")
    }

    if *DEBUG {
        fmt.Println("debug tests enabled")
    } else {
        fmt.Println("debug tests disabled")
    }

    var k, m uint8
    if *k_arg > 255 || *m_arg > 255 {
        log.Fatal("k and m must be >= 255")
    } else {
        k = uint8(*k_arg)
        m = uint8(*m_arg)
        fmt.Println("k =", k)
        fmt.Println("m =", m)
    }

    repeatGenome := GenerateRepeatGenome(genomeName, k, m)

    if *DEBUG {

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

        //fmt.Println("classTree.getLCA(classTree.ClassNodes['DNA/TcMar-Mariner'], classTree.ClassNodes['DNA/TcMar-Tc1']):", classTree.getLCA(classTree.ClassNodes["DNA/TcMar-Mariner"], classTree.ClassNodes["DNA/TcMar-Tc1"]))

        fmt.Println()
        fmt.Println("min(5, 7):", min(5, 7))
        fmt.Println("max64(int64(5), int64(7)):", max64(int64(5), int64(7)))
        fmt.Println()

        offset, thisMin := getMinimizer("tgctcctgtcatgcatacgcaggtcatgcat", 15)
        fmt.Println("getMinimizer('tgctcctgtcatgcatacgcaggtcatgcat', 15) offset :", offset)
        printSeqInt(thisMin, 15)
        fmt.Println()
    }
}
