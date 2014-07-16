package main

import (
    "bytes"
    "flag"
    "fmt"
    "io/ioutil"
    "jh/repeatgenome"
    "os"
    "runtime/pprof"
)

func checkError(err error) {
    if err != nil {
        panic(err)
    }
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

func main() {

    if len(os.Args) < 2 {
        fmt.Println("arg error - usage: ./minimize <flags> <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[len(os.Args)-1]

    cpuProfile := flag.Bool("cpuprof", false, "write cpu profile to file <genomeName>.cpuprof")
    memProfile := flag.Bool("memprof", false, "write memory profile to <genomeName>.memprof")
    debug := flag.Bool("debug", false, "run and print debugging tests")
    noKraken := flag.Bool("no_kraken", false, "don't generate Kraken data structure")
    dontWriteMins := flag.Bool("no_write_kraken", false, "don't write the Kraken data to file")
    writeJSON := flag.Bool("json", false, "write JSON representation of class tree to <genomeName>.classtree.json")
    k_arg := flag.Uint("k", 31, "kmer length")
    m_arg := flag.Uint("m", 15, "minimizer length")
    flag.Parse()

    if *cpuProfile {
        fmt.Println("CPU profiler enabled")
        fmt.Println()
        os.Mkdir("profiles", os.ModeDir)
        f, err := os.Create("profiles/" + genomeName + ".cpuprof")
        if err != nil {
            panic(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    }

    if *memProfile {
        fmt.Println("memory profiler enabled")
        fmt.Println()
    }

    if *debug {
        fmt.Println("debug tests enabled")
        fmt.Println()
    }

    genKraken := !(*noKraken)
    writeKraken := genKraken && !*dontWriteMins
    if !genKraken {
        fmt.Println("Kraken data structure disabled")
        fmt.Println()
    } else if !writeKraken {
        fmt.Println("Writing Kraken data to file disabled")
        fmt.Println()
    }


    if *writeJSON && genKraken {
        fmt.Println("class tree writeJSON write enabled")
        fmt.Println()
    } else if *writeJSON && !genKraken {
        fmt.Println("class tree writeJSON write enabled")
        fmt.Println("WARNING: you are writing the class tree writeJSON without minimizing - all node sizes will be 0")
        fmt.Println()
    }

    var k, m uint8
    if *k_arg > 255 || *m_arg > 255 {
        panic("k and m must be >= 255")
    } else {
        k = uint8(*k_arg)
        m = uint8(*m_arg)
        fmt.Println("k =", k)
        fmt.Println("m =", m)
        fmt.Println()
    }

    rgFlags := repeatgenome.Flags{*debug, *cpuProfile, *memProfile, genKraken, writeKraken, *writeJSON}
    repeatGenome := repeatgenome.Generate(genomeName, k, m, rgFlags)

    workingDirName, err := os.Getwd()
    checkError(err)
    readsDirName := workingDirName + "/" + genomeName + "-reads"
    currDir, err := os.Open(readsDirName)
    checkError(err)
    fileinfos, err := currDir.Readdir(-1)
    checkError(err)
    processedFiles := []os.FileInfo{}
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 5 && fileinfo.Name()[len(fileinfo.Name())-5 : ] == ".proc" {
            processedFiles = append(processedFiles, fileinfo)
        }
    }
    readsBytes := [][]byte{}
    for _, fileinfo := range processedFiles {
        _, theseReadsBytes := fileLines(readsDirName + "/" + fileinfo.Name())
        for _, bytesLine := range theseReadsBytes {
            readsBytes = append(readsBytes, bytesLine)
        }
    }

    readChan := make(chan string)

    go func() {
        for _, readBytes := range readsBytes {
            readChan <- string(readBytes)
        }
        close(readChan)
    }()

    repeatGenome.ClassifyReads(readChan)
    fmt.Println(repeatGenome.Name, "successfully generated - exiting")
}
