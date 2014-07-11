package repeatgenome

/*
   A collection of trivial functions used in other source files of the repeatgenome package.
*/

import (
    "bytes"
    "io/ioutil"
    "log"
    //"reflect"
    "sync"
)

var err error

func checkError(err error) {
    if err != nil {
        log.Fatal(err)
    }
}

func merge(cs [](chan interface{})) <-chan interface{} {
    var wg sync.WaitGroup
    //elemType := reflect.TypeOf(cs).Elem()
    //chanType := reflect.ChanOf(RecvDir, elemType)
    out := make(chan interface{})

    // Start an output goroutine for each input channel in cs.  output
    // copies values from c to out until c is closed, then calls wg.Done.
    wg.Add(len(cs))
    for _, c := range cs {
        go func(c chan interface{}) {
            for n := range c {
                    out <- n
            }
            wg.Done()
        }(c)
    }

    // Start a goroutine to close out once all the output goroutines are
    // done.  This must start after the wg.Add call.
    go func() {
        wg.Wait()
        close(out)
    }()
    return out
}

// returns the number of lines and a slice of the lines
func lines(byteSlice []byte) [][]byte {
    var lines [][]byte = bytes.Split(byteSlice, []byte{'\n'})
    // drop the trailing newlines
    newline := []byte("\n")
    for lastLine := lines[len(lines)-1]; len(lines) > 0 && (len(lastLine) == 0 || bytes.Equal(lastLine, newline)); lastLine = lines[len(lines)-1] {
        lines = lines[:len(lines)-1]
    }
    return lines
}

func fileLines(filepath string) (err error, linesBytes [][]byte) {
    rawBytes, err := ioutil.ReadFile(filepath)
    if err != nil {
        return err, nil
    } else {
        return nil, lines(rawBytes)
    }
}

/*
func readSimSeqReads(filepath string) (err error, numReads uint64, simReads []Seq) {
    err, 
}
*/

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

func minU64(a, b uint64) uint64 {
    if a < b {
        return a
    } else {
        return b
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
    return bytes.Compare(kmers[i][:8], kmers[j][:8]) == -1
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
