main: repeatgenome.go helper-rg.go stats-rg.go run-rg.go io-rg.go seq-rg.go
	go install jh/repeatgenome
	go build run-rg.go

gcc: repeatgenome.go helper-rg.go stats-rg.go run-rg.go io-rg.go seq-rg.go
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o jh/repeatgenome.o helper-rg.go io-rg.go stats-rg.go seq-rg.go repeatgenome.go
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o run-rg.o run-rg.go
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -o run-rg-gcc run-rg.o jh/repeatgenome.o
	rm jh/repeatgenome.o run-rg.o
