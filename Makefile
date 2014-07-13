main: repeatgenome.go helper-rg.go stats.go run-rg.go
	go install jh/repeatgenome
	go build run-rg.go

gcc: repeatgenome.go helper-rg.go stats.go run-rg.go
	gccgo -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o jh/repeatgenome.o helper-rg.go stats.go repeatgenome.go
	gccgo -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o run-rg.o run-rg.go
	gccgo -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -o run-rg-gcc run-rg.o jh/repeatgenome.o
