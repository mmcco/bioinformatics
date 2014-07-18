package repeatgenome

import (
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
