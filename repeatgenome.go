package repeatgenome

/*
   WARNING!!! This program is currently under development and may be buggy or broken.

   A barebones (at the moment) Go script for parsing and minimizing RepeatMasker output files alongside FASTA reference genomes.

   This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
       * a RepeatMasker library containing:
           - the match library (e.g. "dm3.fa.out")
           - the alignment information (e.g. "dm3.fa.align")
       * one or more reference genome files in FASTA format with the suffix ".fa"

   Premature commenting is the root of all evil, and I have sinned. Please read comments skeptically - they will eventually be audited.

   For portability's sake, the flags should be used as args to Generate() rather than globals.

   Minimizers could be stored as uint32's.

   We should test a version that doesn't cache minimizers, as that seems to be a needless bottleneck. It could also be conditional on the number of CPUs available.

   Minimizers are currently not written in any order. This is for memory efficiency, and can be changed in the future.

   To work around the missing repeat ID RepeatMasker bug, we should probably do a lookup on the repeat ID in a map before failing.

   I should probably change some variable names, like repeatGenome, to less verbose variants, and use more aliases.

   Ranges should be changed to use actual values instead of indexes.

   Slice sizes should be specified in the make() call when the size is known.

   seqToInt and revCompToInt need minor cleanup and a potential name-change.

   Should reconsider what is a pointer and what is directly referenced

   All sequences containing Ns are currently ignored.

   We should consider taking end minimizers once the code base is more mature.

   We should also review how to deal with m <= len(match) < k.

   For caching efficiency, we should change the minimizer data structure to a map-indexed 1D slice of Kmers (not *Kmers). (This technique originated in Kraken.)

   Int sizes should be reviewed for memory efficiency.

   The sole command line argument is the name of the reference genome (e.g. "dm3").
*/

import (
    "bytes"
    "fmt"
    "io/ioutil"
    "log"
    //"mapset"
    "os"
    "runtime"
    "runtime/pprof"
    "sort"
    "strconv"
    "strings"
    "sync"
    "unsafe"
)

type ParseFlags struct {
    Debug       bool
    CPUProfile  bool
    MemProfile  bool
    Minimize    bool
    WriteKraken bool
    WriteJSON   bool
}

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
// Match.InsertionID - a numerical ID for the repeat type (starts at 1)
//     Match.Fullname - not in parsed data file - simply repeatClass concatenated - used for map indexing
type Match struct {
    SW_Score    int32
    PercDiv     float64
    PercDel     float64
    PercIns     float64
    SeqName     string
    SeqStart    uint64
    SeqEnd      uint64
    SeqRemains  uint64
    IsRevComp   bool
    RepeatClass []string
    // in weird cases, RepeatStart can be negative, so they must be signed
    RepeatStart   int64
    RepeatEnd     int64
    RepeatRemains int64
    InsertionID   uint64

    // these are generated, not parsed
    RepeatName string
    ClassNode  *ClassNode
    Repeat     *Repeat
}

type RepeatGenome struct {
    Name       string
    ParseFlags ParseFlags
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    chroms       map[string](map[string]string)
    K            uint8
    M            uint8
    Kmers        []Kmer
    OffsetsToMin map[uint64]uint64
    SortedMins   Uint64Slice
    Matches      Matches
    ClassTree    ClassTree
    Repeats      Repeats
    RepeatMap    map[string]uint64
}

type ClassTree struct {
    // maps all class names to a pointer to their node struct
    // we must use pointers because of this foible in golang: https://code.google.com/p/go/issues/detail?id=3117
    // if we didn't use pointers, we'd have to completely reassign the struct when adding parents, children, etc.
    ClassNodes map[string](*ClassNode)
    NodesByID  []*ClassNode
    // a pointer to the the class tree's root, used for recursive descent etc.
    // we explicitly create the root (because RepeatMatcher doesn't)
    Root *ClassNode
}

// can store a kmer where k <= 32
// the value of k is not stored in the struct, but rather in the RepeatGenome, for memory efficiency
// first eight bits are the int representation of the sequence
// the last two are the LCA ID
type Kmer [10]byte

// can store a sequence of length <= 65,535
type Seq struct {
    Bases  []byte
    Length uint16
}

type Repeat struct {
    // assigned in simple incremented order starting from 1
    // they are therefore not compatible across genomes
    // we give root ID = 0
    ID uint64
    // a list containing the repeat's ancestry path, from top down
    // root is implicit, and is therefore excluded from the list
    ClassList []string
    ClassNode *ClassNode
    Name  string
    NumInstances uint64
}

type ClassNode struct {
    Name     string
    ID       uint16
    Class    []string
    Parent   *ClassNode
    Children []*ClassNode
    IsRoot   bool
}

// type synonyms, necessary to implement interfaces (e.g. sort) and methods
type Kmers []*Kmer
type MinMap map[uint64]Kmers
type Repeats []Repeat
type Matches []Match

//type Chroms map[string](map[string][]byte)

type MinPair struct {
    KmerInt   Kmer
    Minimizer uint64
}

type ThreadResponse struct {
    KmerInt   uint64
    Minimizer uint64
    Relative *ClassNode
}

type MinCache struct {
    sync.Mutex
    Cache map[uint64]uint64
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
        switch seq[len(seq)-(i+1)] {
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

func parseMatches(genomeName string) Matches {
    // "my_genome_name"  ->  "my_genome_name/my_genome_name.fa.out"
    filepath := strings.Join([]string{genomeName, "/", genomeName, ".fa.out"}, "")
    err, matchLines := fileLines(filepath)
    checkError(err)
    // drop header
    matchLines = matchLines[3:]

    var matches Matches
    var sw_Score int64

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
        match.InsertionID, err = strconv.ParseUint(rawVals[14], 10, 64)
        checkError(err)

        // necessary swaps to convert reverse complement repeat indexes to positive-strand indexes
        if match.IsRevComp {
            match.RepeatStart = match.RepeatRemains
            match.RepeatEnd = match.RepeatStart
            match.RepeatRemains = match.RepeatRemains + (match.RepeatEnd - match.RepeatStart)
        }

        // decrement match.SeqStart and match.RepeatStart so that they work from a start index of 0 rather than 1
        // that way, we can use them without modification in slices
        match.SeqStart--
        match.RepeatStart--

        // "Other" and "Unknown" classes are heirarchically meaningless and really just mean "root", so we remove them
        if match.RepeatClass[0] == "Other" || match.RepeatClass[0] == "Unknown" {
            match.RepeatClass = match.RepeatClass[1:]
        }

        match.RepeatName = strings.Join(match.RepeatClass, "/")

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
            err, seqLines := fileLines(chromFilepath)
            checkError(err)

            // maps each sequence name in this chrom to a slice of its sequence's lines
            // the list is concatenated at the end for efficiency's sake
            seqMap := make(map[string][][]byte)
            var numLines uint64 = uint64(len(seqLines))
            var seqName string
            var i uint64
            for i = 0; i < numLines; i++ {
                seqLine := bytes.TrimSpace(seqLines[i])
                if seqLine[0] == byte('>') {
                    seqName = string(bytes.TrimSpace(seqLine[1:]))
                    if !warned && seqName != chromFilename[:len(chromFilename)-3] {
                        fmt.Println("WARNING: reference genome is two-dimensional, containing sequences not named after their chromosome.")
                        fmt.Println("Because RepeatMasker supplied only one-dimensional indexing, this may cause unexpected behavior or program failure.")
                        fmt.Println("seqName:", seqName, "\tlen(seqName):", len(seqName))
                        fmt.Println("chrom name:", chromFilename[:len(chromFilename)-3], "\tlen(chrom name):", len(chromFilename)-3)
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

func (repeatGenome *RepeatGenome) Size() uint64 {
    var size uint64 = 0

    for _, seqs := range repeatGenome.chroms {
        for _, seq := range seqs {
            size += uint64(len(seq))
        }
    }
    return size
}

func Generate(genomeName string, k, m uint8, parseFlags ParseFlags) *RepeatGenome {
    // we popoulate the RepeatGenome mostly with helper functions
    // we should consider whether it makes more sense for them to alter the object directly, than to return their results
    repeatGenome := new(RepeatGenome)
    repeatGenome.Name = genomeName
    repeatGenome.ParseFlags = parseFlags
    repeatGenome.chroms = parseGenome(genomeName)
    repeatGenome.Matches = parseMatches(genomeName)
    repeatGenome.getRepeats()
    repeatGenome.getClassTree()
    repeatGenome.K = k
    repeatGenome.M = m

    if repeatGenome.ParseFlags.Minimize {
        // calling the parallel minimizer and writing the result
        repeatGenome.getKrakenSlice()
    }

    if repeatGenome.ParseFlags.WriteJSON {
        repeatGenome.WriteClassJSON(false, false)
    }

    if repeatGenome.ParseFlags.Debug {

        fmt.Println()
        for k, v := range repeatGenome.chroms {
            for k_, v_ := range v {
                fmt.Printf("chrom: %s\tseq: %s\t%s...%s\n", k, k_, v_[:20], v_[len(v_)-20:])
            }
        }
        fmt.Println()

        fmt.Println()
        fmt.Println("number of chromosomes parsed:", len(repeatGenome.chroms))
        fmt.Println()

        fmt.Println("total number of bases in genome:", repeatGenome.Size())

        repeatGenome.ClassTree.PrintBranches()
        fmt.Println()
        fmt.Println("number of ClassNodes:", len(repeatGenome.ClassTree.ClassNodes))
        fmt.Println()

        fmt.Println("repeatGenome.ClassTree.getLCA(repeatGenome.ClassTree.ClassNodes['DNA/TcMar-Mariner'], repeatGenome.ClassTree.ClassNodes['DNA/TcMar-Tc1']):", repeatGenome.ClassTree.getLCA(repeatGenome.ClassTree.ClassNodes["DNA/TcMar-Mariner"], repeatGenome.ClassTree.ClassNodes["DNA/TcMar-Tc1"]))

        fmt.Println()
        fmt.Println("min(5, 7):", min(5, 7))
        fmt.Println("max64(int64(5), int64(7)):", max64(int64(5), int64(7)))
        fmt.Println()

        testSeq := "atgtttgtgtttttcataaagacgaaagatg"
        offset, thisMin := getMinimizer(seqToInt(testSeq), uint8(len(testSeq)), 15)
        fmt.Println("getMinimizer('tgctcctgtcatgcatacgcaggtcatgcat', 15) offset :", offset)
        printSeqInt(thisMin, 15)
        fmt.Println()

        fmt.Printf("Kmer struct size: %d\n", unsafe.Sizeof(Kmer{}))
    }

    return repeatGenome
}

func (repeatGenome *RepeatGenome) getRepeats() {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID
    repeatGenome.Repeats = Repeats{}

    // maps a repeat's category to its ID
    repeatGenome.RepeatMap = make(map[string]uint64)

    for _, match := range repeatGenome.Matches {
        // don't bother overwriting
        repeatID, exists := repeatGenome.RepeatMap[match.RepeatName]
        if exists {
            repeatGenome.Repeats[repeatID].NumInstances++
            match.Repeat = &repeatGenome.Repeats[repeatID]
        } else {
            repeat := Repeat{}
            repeat.ID = uint64(len(repeatGenome.Repeats))
            repeat.ClassList = match.RepeatClass
            repeat.Name = match.RepeatName
            repeat.NumInstances = 1

            repeatGenome.Repeats = append(repeatGenome.Repeats, repeat)
            repeatGenome.RepeatMap[repeat.Name] = repeat.ID

            match.Repeat = &repeat
        }
    }
}

func getMinimizer(kmer uint64, k, m uint8) (uint8, uint64) {
    if m > k || m < 1 {
        panic("getMinimizer(): m must be <= k and > 0")
    }

    // stores the index of the leftmost base included in the minimizer
    var currOffset uint8 = 0
    numExtraBases := 32 - k
    revCompKmer := intRevComp(kmer, k)
    currMin := kmer >> uint64(64-2*(32-k+m))
    possMin := currMin
    numHangingBases := k - m
    var i uint8
    // i is the index of the offset
    for i = 0; i <= numHangingBases; i++ {
        // overflow off the first excluded base
        possMin = kmer << (2 * (numExtraBases + i))
        // return to proper alignment
        possMin >>= 64 - 2*m

        if possMin < currMin {
            currMin = possMin
            currOffset = i
        }

        possMin = revCompKmer << (2 * (numExtraBases + i))
        possMin >>= 64 - 2*m

        if possMin < currMin {
            currMin = possMin
            currOffset = numHangingBases - i
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
        default:
            panic("byte other than 'a', 'c', 'g', or 't' supplied to revComp")
        }
    }
    return string(revCompSeq)
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

func (repeatGenome *RepeatGenome) getClassTree() {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    tree := &repeatGenome.ClassTree
    tree.ClassNodes = make(map[string](*ClassNode))
    tree.NodesByID = make([]*ClassNode, 0)
    // would be prettier if expanded
    tree.Root = &ClassNode{"root", 0, []string{"root"}, nil, nil, true}
    tree.ClassNodes["root"] = tree.Root
    tree.NodesByID = append(tree.NodesByID, tree.Root)

    for _, repeat := range repeatGenome.Repeats {
        // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
        for j := 1; j <= len(repeat.ClassList); j++ {
            thisClass := repeat.ClassList[:j]
            thisClassName := strings.Join(thisClass, "/")
            _, keyExists := tree.ClassNodes[thisClassName]
            if !keyExists {
                if len(tree.NodesByID) > 65534 {
                    panic("RepeatGenome.getClassTree(): more than 65,536 class nodes - ID is overflowed")
                }
                classNode := new(ClassNode)
                classNode.Name = thisClassName
                classNode.ID = uint16(len(tree.NodesByID))
                classNode.Class = thisClass
                classNode.IsRoot = false

                tree.ClassNodes[thisClassName] = classNode
                tree.NodesByID = append(tree.NodesByID, classNode)
                // first case handles primary classes, as root is implicit and not listed in thisClass
                if j == 1 {
                    classNode.Parent = tree.Root
                } else {
                    classNode.Parent = tree.ClassNodes[strings.Join(thisClass[:len(thisClass)-1], "/")]
                }
                if classNode.Parent.Children == nil {
                    classNode.Parent.Children = make([]*ClassNode, 0)
                }
                classNode.Parent.Children = append(classNode.Parent.Children, tree.ClassNodes[thisClassName])
            }
        }
    }

    // MUST NOT USE RANGE - the struct will be copied!
    for i := 0; i < len(repeatGenome.Repeats); i++ {
        repeat := &repeatGenome.Repeats[i]
        repeat.ClassNode = tree.ClassNodes[repeat.Name]

        if repeat.ClassNode == nil {
            fmt.Println(repeat.Name)
            log.Fatal("getClassTree(): nil Repeat.ClassNode")
        }
    }

    // MUST NOT USE RANGE - the struct will be copied!
    for i := 0; i < len(repeatGenome.Matches); i++ {
        match := &repeatGenome.Matches[i]
        match.ClassNode = tree.ClassNodes[match.RepeatName]

        if match.ClassNode == nil {
            fmt.Println(match.RepeatName)
            log.Fatal("getClassTree(): nil Match.ClassNode")
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

// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func (repeatGenome *RepeatGenome) minimizeThread(matchStart, matchEnd uint64, c chan ThreadResponse) {
    var match *Match
    k := repeatGenome.K
    k_ := uint64(k)
    m := repeatGenome.M
    //m_ := uint64(m)
    var seq string
    var x, matchLen, kmerInt, thisMin uint64

    for i := matchStart; i < matchEnd; i++ {
        match = &repeatGenome.Matches[i]
        seq = repeatGenome.chroms[match.SeqName][match.SeqName]
        matchLen = match.SeqEnd - match.SeqStart
        // for now, we will ignore matches too short to be traditionally minimized
        if matchLen < k_ {
            continue
        }

    KmerLoop:
        for j := match.SeqStart; j <= match.SeqEnd-k_; j++ {
            // we begin by skipping any kmers containing n's
            // we start checking from the end for maximum skipping efficiency
            for x = j + k_ - 1; x >= j; x-- {
                if seq[x] == byte('n') {
                    j += x
                    continue KmerLoop
                }
                // prevents negative overflow - a bit of a hack, but seqs can be big, so we need uint64's capacity
                if x == 0 {
                    break
                }
            }

            kmerInt = seqToInt(seq[j : j+k_])
            _, thisMin = getMinimizer(kmerInt, k, m)
            if match.ClassNode == nil {
                fmt.Println("current match's repeat:", match.RepeatName)
                panic("minimizeThread(): match has nil ClassNode")
            }
            c <- ThreadResponse{kmerInt, thisMin, match.ClassNode}
        }
    }

    close(c)
}

func (repeatGenome *RepeatGenome) getKrakenSlice() {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    numCPU := runtime.NumCPU()
    if repeatGenome.ParseFlags.Debug {
        fmt.Printf("using %d CPUs\n", numCPU)
    }
    runtime.GOMAXPROCS(numCPU)
    var mStart, mEnd uint64

    numKmers := repeatGenome.numKmers()
    fmt.Printf("expecting >= %d million kmers\n", numKmers/1000000)

    var threadChans [](chan ThreadResponse)
    for i := 0; i < numCPU; i++ {
        mStart = uint64(i * len(repeatGenome.Matches) / numCPU)
        mEnd = uint64((i + 1) * len(repeatGenome.Matches) / numCPU)
        c := make(chan ThreadResponse, 1000)
        threadChans = append(threadChans, c)
        go repeatGenome.minimizeThread(mStart, mEnd, c)
    }

    // below is the atomic section of minimizing, which is parallel
    // this seems to be the rate-limiting section, as threads use only ~6-7 CPU-equivalents
    // it should therefore be optimized before other sections
    var kmer *Kmer
    var kmersProcessed, kmerInt, minimizer uint64
    var lca_ID uint16
    var lca, relative *ClassNode
    var exists bool
    minCache := make(map[uint64]uint64)
    kmerMap := make(map[uint64]*Kmer)

    for response := range merge(threadChans) {
        if kmersProcessed%5000000 == 0 {
            fmt.Println(kmersProcessed/1000000, "million kmers processed")
        }
        kmersProcessed++

        kmerInt, minimizer, relative = response.KmerInt, response.Minimizer, response.Relative
        minCache[kmerInt] = minimizer

        kmer, exists = kmerMap[kmerInt]
        if exists {
            lca_ID = *(*uint16)(unsafe.Pointer(&kmer[8]))
            lca = repeatGenome.ClassTree.getLCA(repeatGenome.ClassTree.NodesByID[lca_ID], relative)
            *(*uint16)(unsafe.Pointer(&kmer[8])) = lca.ID
        } else {
            kmer = new(Kmer)
            *(*uint64)(unsafe.Pointer(&kmer[0])) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = relative.ID
            kmerMap[kmerInt] = kmer
        }
    }

    fmt.Println("all minimizers generated")
    fmt.Println(len(kmerMap), "unique kmers generated")

    var kmers Kmers
    var thisMin uint64
    minMap := make(map[uint64]Kmers)

    for kmerInt, kmer := range kmerMap {
        minimizer = minCache[kmerInt]
        kmers, exists = minMap[minimizer]
        if exists {
            minMap[minimizer] = append(kmers, kmer)
        } else {
            minMap[minimizer] = Kmers{kmer}
        }
        delete(kmerMap, kmerInt)
    }

    // manual deletion to save memory
    kmerMap = nil
    runtime.GC()

    if repeatGenome.ParseFlags.WriteKraken {
        err = repeatGenome.WriteMins(minMap)
        checkError(err)
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
        thisMin = repeatGenome.SortedMins[i]
        repeatGenome.OffsetsToMin[thisMin] = currOffset
        currOffset += uint64(len(minMap[thisMin]))
        for j := range minMap[thisMin] {
            repeatGenome.Kmers = append(repeatGenome.Kmers, *minMap[thisMin][j])
        }
        delete(minMap, thisMin)
    }

    if repeatGenome.ParseFlags.MemProfile {
        os.Mkdir("profiles", os.ModeDir)
        f, err := os.Create("profiles/" + repeatGenome.Name + ".memprof")
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

    splitOnN := func(c rune) bool { return c == 'n' }

    for i := range repeatGenome.Matches {
        match = &repeatGenome.Matches[i]
        seq = repeatGenome.chroms[match.SeqName][match.SeqName][match.SeqStart:match.SeqEnd]
        seqs = strings.FieldsFunc(seq, splitOnN)
        for j := range seqs {
            if len(seqs[j]) >= k {
                numKmers += uint64(len(seqs[j]) - k + 1)
            }
        }
    }
    return numKmers
}
