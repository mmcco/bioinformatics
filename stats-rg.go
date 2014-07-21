package repeatgenome

import (
    "fmt"
    "unsafe"
)

func (repeatGenome *RepeatGenome) Size() uint64 {
    var numBases uint64 = 0
    for _, seqs := range repeatGenome.chroms {
        for _, seq := range seqs {
            numBases += uint64(len(seq))
        }
    }
    return numBases
}

func (rg *RepeatGenome) KmersGBSize() float64 {
    return (float64(len(rg.Kmers))/1000000000) * float64(unsafe.Sizeof(Kmer{}))
}

func (rg *RepeatGenome) PercentRepeats() float64 {
    var totalBases uint64 = rg.Size()

    var repeatBases uint64 = 0
    for _, match := range rg.Matches {
        repeatBases += match.SeqEnd - match.SeqStart
    }

    return 100 * (float64(repeatBases) / float64(totalBases))
}

// written for the PercentTrueClassification() below
// determines whether a read overlaps any repeat instances in the given ClassNode's subtree
func recNodeSearch(classNode *ClassNode, readSAM ReadSAM) bool {
    if classNode != nil && classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            // must compute where the read ends
            endInd := readSAM.StartInd + uint64(len(readSAM.Seq))
            if readSAM.SeqName == match.SeqName && readSAM.StartInd < match.SeqEnd && endInd > match.SeqStart {
                return true
                // below logic is for checking for at least rg.K overlap
                /*
                var overlap uint64 := readSAM.SeqEnd - match.SeqStart
                if readSAM.SeqStart > match.SeqStart {
                    overlap -= readSAM.SeqStart - match.SeqStart
                }
                if overlap >= uint64(rg.K) {
                    return true
                }
                */
            }
        }
    }
    if classNode != nil && classNode.Children != nil {
        for _, child := range classNode.Children {
            if recNodeSearch(child, readSAM) {
                return true
            }
        }
    }
    return false
}

func TestNodeSearch(classNode *ClassNode, readSAM ReadSAM) bool {
    if classNode == nil || (classNode.Name != "root" && classNode.Name != "Satellite" && classNode.Name != "Satellite/HETRP_DM") {
        return false
    }
    fmt.Println("testing", classNode.Name)
    if classNode != nil && classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            // must compute where the read ends
            endInd := readSAM.StartInd + uint64(len(readSAM.Seq))
            fmt.Printf("testing match %s[%d:%d] against %s[%d:%d]\n", match.SeqName, match.SeqStart, match.SeqEnd, readSAM.SeqName, readSAM.StartInd, endInd)
            if readSAM.SeqName == match.SeqName && readSAM.StartInd < match.SeqEnd && endInd > match.SeqStart {
                fmt.Println("true")
                return true
                // below logic is for checking for at least rg.K overlap
                /*
                var overlap uint64 := readSAM.SeqEnd - match.SeqStart
                if readSAM.SeqStart > match.SeqStart {
                    overlap -= readSAM.SeqStart - match.SeqStart
                }
                if overlap >= uint64(rg.K) {
                    return true
                }
                */
            }
        }
    } else {
        fmt.Println("no classNode or no classNode.Repeat")
    }
    if classNode != nil && classNode.Children != nil {
        for _, child := range classNode.Children {
            if TestNodeSearch(child, readSAM) {
                fmt.Println("true")
                return true
            } else {
                fmt.Println("(child false)")
            }
        }
    } else {
        fmt.Println("no classNode or no classNode.Children")
    }
    fmt.Println("false")
    return false
}

// we currently use the simple metric that the read and one of the repeat's instances overlap at all
func (rg *RepeatGenome) PercentTrueClassifications(responses []ReadSAMResponse) float64 {
    var classifications, correctClassifications uint64 = 0, 0

    for _, resp := range responses {
        if resp.ClassNode != nil {
            classifications++
        }
        if recNodeSearch(resp.ClassNode, resp.ReadSAM) {
            correctClassifications++
        }
    }
    return 100 * (float64(correctClassifications) / float64(classifications))
}