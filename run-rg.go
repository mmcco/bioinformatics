package main

import (
    "flag"
    "fmt"
    "jh/repeatgenome"
    "os"
    "runtime/pprof"
)

func main() {

    if len(os.Args) < 2 {
        fmt.Println("arg error - usage: ./minimize <flags> <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[len(os.Args)-1]

    cpuProfile := flag.Bool("cpuprof", false, "write cpu profile to file <genomeName>.cpuprof")
    memProfile := flag.Bool("memprof", false, "write memory profile to <genomeName>.memprof")
    debug := flag.Bool("debug", false, "run and print debugging tests")
    dont_minimize := flag.Bool("nomin", false, "don't generate Kraken data structure")
    writeJSON := flag.Bool("json", false, "write writeJSON representation of class tree to <genomeName>.classtree.json")
    k_arg := flag.Uint("k", 31, "kmer length")
    m_arg := flag.Uint("m", 15, "minimizer length")
    flag.Parse()

    if *cpuProfile {
        fmt.Println("CPU profiler enabled")
        fmt.Println()
        f, err := os.Create(genomeName + ".cpuprof")
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

    var minimize bool
    if *dont_minimize {
        minimize = false
        fmt.Println("Kraken data structure disabled")
        fmt.Println()
    } else {
        minimize = true
    }

    if *writeJSON && minimize {
        fmt.Println("class tree writeJSON write enabled")
        fmt.Println()
    } else if *writeJSON && !minimize {
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

    parseFlags := repeatgenome.ParseFlags{*debug, *cpuProfile, *memProfile, minimize, *writeJSON}
    repeatGenome := repeatgenome.Generate(genomeName, k, m, parseFlags)
    fmt.Println(repeatGenome.Name, "successfully generated - exiting")
}
