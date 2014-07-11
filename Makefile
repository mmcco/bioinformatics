main: repeatgenome.go helper-rg.go stats.go run-rg.go
	go install jh/repeatgenome
	go build run-rg.go
