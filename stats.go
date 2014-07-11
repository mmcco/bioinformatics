package repeatgenome

func (repeatGenome *RepeatGenome) size() uint64 {
    var numBases uint64 = 0
    for _, seqs := range repeatGenome.chroms {
        for _, seq := range seqs {
            numBases += uint64(len(seq))
        }
    }
    return numBases
}
