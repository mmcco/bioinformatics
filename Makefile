package_dir = $(GOPATH)/src/jh/repeatgenome/

main: $(package_dir)repeatgenome.go $(package_dir)helper-rg.go $(package_dir)stats-rg.go $(package_dir)io-rg.go $(package_dir)seq-rg.go run-rg.go 
	go install jh/repeatgenome
	go build run-rg.go

gcc: $(package_dir)repeatgenome.go $(package_dir)helper-rg.go $(package_dir)stats-rg.go $(package_dir)io-rg.go $(package_dir)seq-rg.go run-rg.go 
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o jh/repeatgenome.o helper-rg.go io-rg.go stats-rg.go seq-rg.go repeatgenome.go
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -c -o run-rg.o run-rg.go
	gccgo -Wall -static -L /scratch0/langmead-fs1/shared/go/pkg/linux_amd64/ -lpthread -O2 -o run-rg-gcc run-rg.o jh/repeatgenome.o
	rm jh/repeatgenome.o run-rg.o
