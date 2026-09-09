.PHONY: test test-all test-free test-wolfram

test: test-all

test-all: test-free test-wolfram

test-free:
	./test/run-free-tests.sh

test-wolfram:
	./test/run-wolfram-tests.sh
